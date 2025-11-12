//! Search algorithms for finding primers in DNA sequences.
//!
//! This module contains implementations of different search algorithms for finding
//! primer sequences within reads. Each search algorithm returns an `Option<PrimerMatch>`
//! to `search_paired_primers` where:
//! - `Some(PrimerMatch)` indicates a match was found with the best match details.
//!   The `PrimerMatch` contains `start` and `end` positions (0-based) that are used
//!   by `search_paired_primers` to determine trim positions for removing primers from reads.
//!   The positions indicate where the primer was found in the read sequence.
//! - `None` indicates no match was found within the specified constraints.
//!
//! The search algorithms return the same type (`PrimerMatch`) and can be used
//! by `search_paired_primers` to locate primers at the start and end of reads.

use clap::ValueEnum;

/// Algorithm selection for primer matching
#[derive(ValueEnum, Clone, Debug, Copy)]
pub enum Algorithm {
    /// Use Myers bit-parallel algorithm for approximate matching
    Myers,
    /// Use Sassy SIMD-accelerated search (fastest, requires AVX2/NEON)
    Sassy,
    /// Use BNDM for exact matching (fastest for short exact matches)
    Bndm,
}

impl Default for Algorithm {
    fn default() -> Self {
        Algorithm::Sassy
    }
}
use crate::preparing_input::build_myers_matcher;
use crate::primer_search::SearchConfig;
use bio::pattern_matching::bndm::BNDM;
use bio::pattern_matching::myers::long::Myers;
use sassy::{Searcher, profiles::Iupac};

/// Result of a primer search
/// Contains the positions where a primer was found in a read sequence
#[derive(Debug, Clone)]
pub struct PrimerMatch {
    pub start: usize,
    pub end: usize,
}

/// Find primer in read using Myers bit-parallel algorithm with IUPAC support
/// If myers_cache is provided, uses the pre-built matcher; otherwise builds one on the fly
fn find_primer_myers(
    cfg: &SearchConfig,
    read: &str,
    primer: &str,
    myers_cache: Option<&mut Myers<u64>>,
) -> Option<PrimerMatch> {
    // Search in the first 'window' bases
    let window = cfg.window.min(read.len());
    let region_bytes = read[..window].as_bytes();
    let min_overlap_len = cfg.min_overlap.min(primer.len());

    // Calculate max distance based on error rate
    // For a primer match, we want matches within error_rate of the alignment length
    let max_dist = ((primer.len() as f64) * cfg.error_rate).ceil() as usize;

    // Use cached matcher if available, otherwise build one temporarily
    let mut owned_matcher;
    let myers: &mut Myers<u64> = if let Some(cached) = myers_cache {
        cached
    } else {
        owned_matcher = build_myers_matcher(primer.as_bytes());
        &mut owned_matcher
    };

    // Find all matches within max distance using find_all to get (start, end, distance)
    let mut best_match: Option<PrimerMatch> = None;
    let mut best_distance = usize::MAX;

    // Iterate through all matches to find the best one
    for (start, end, distance) in myers.find_all(region_bytes, max_dist) {
        let overlap_len = end - start;

        // Check minimum overlap requirement
        if overlap_len < min_overlap_len {
            continue;
        }

        // Calculate actual error rate for this match
        let error_rate = (distance as f64) / (overlap_len as f64);

        // Check error rate threshold
        if error_rate > cfg.error_rate {
            continue;
        }

        if distance >= best_distance {
            continue;
        }

        // This is a valid match!
        // Note: find_all returns end position as exclusive (Rust range convention)
        // We need to adjust to be inclusive for PrimerMatch
        best_match = Some(PrimerMatch {
            start,
            end: end - 1,
        });
        best_distance = distance;
    }

    best_match
}

/// Find primer in read using BNDM for exact matching
/// Searches all concrete variants (expanded from degenerate codes)
fn find_primer_bndm(
    cfg: &SearchConfig,
    read: &str,
    primer_variants: &[String],
) -> Option<PrimerMatch> {
    // Search in the first 'window' bases
    let window = cfg.window.min(read.len());
    let region_bytes = read[..window].as_bytes();
    let min_overlap_len = cfg.min_overlap.min(if primer_variants.is_empty() {
        0
    } else {
        primer_variants[0].len()
    });

    // Try each primer variant (from degenerate expansion)
    for primer_variant in primer_variants {
        let primer_bytes = primer_variant.as_bytes();
        let matcher = BNDM::new(primer_bytes);

        // Find first exact match in the window
        if let Some(pos) = matcher.find_all(region_bytes).next() {
            let match_end = pos + primer_bytes.len();
            let overlap_len = match_end - pos;

            // Check minimum overlap requirement
            if overlap_len >= min_overlap_len {
                return Some(PrimerMatch {
                    start: pos,
                    end: match_end - 1, // Make end inclusive to match Myers convention
                });
            }
        }
    }

    None
}

/// Find primer in read using Sassy SIMD-accelerated search
fn find_primer_sassy(cfg: &SearchConfig, read: &str, primer: &str) -> Option<PrimerMatch> {
    // Search in the first 'window' bases
    let window = cfg.window.min(read.len());
    let region_bytes = &read[..window].as_bytes();

    // Use overhang cost to handle partial primer matches at read boundaries
    // Overhang cost of 0.5 means each overhanging character is penalized by 0.5
    // instead of the full cost - important for truncated primers at contig/read ends
    let overhang_cost = 0.5;
    let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(overhang_cost);
    // Convert error_rate to max allowed edits (generous upper bound)
    let max_edits =
        ((primer.len() as f64) * cfg.error_rate / (1.0 - cfg.error_rate)).ceil() as usize;
    let matches = searcher.search(primer.as_bytes(), region_bytes, max_edits);

    // Find best match (lowest cost)
    let best_match = matches.iter().min_by_key(|m| m.cost).and_then(|m| {
        let overlap_len = m.text_end - m.text_start;
        let min_overlap_len = cfg.min_overlap.min(primer.len());
        let error_rate = (m.cost as f64) / (overlap_len as f64);

        // Check minimum overlap requirement and error rate
        if overlap_len >= min_overlap_len && error_rate <= cfg.error_rate {
            Some(PrimerMatch {
                start: m.text_start,
                end: m.text_end,
            })
        } else {
            None
        }
    });

    best_match
}

/// Find primer in read using degenerate-aware algorithms (Myers, Sassy)
/// Optionally uses a pre-built Myers matcher for the Myers algorithm
pub fn find_primer_degenerate(
    cfg: &SearchConfig,
    read: &str,
    primer: &str,
    myers_cache: Option<&mut Myers<u64>>,
) -> Option<PrimerMatch> {
    match cfg.algorithm {
        Algorithm::Myers => find_primer_myers(cfg, read, primer, myers_cache),
        Algorithm::Sassy => find_primer_sassy(cfg, read, primer),
        Algorithm::Bndm => unreachable!("Exact-match algorithms require pre-expanded variants"),
    }
}

/// Find primer in read using exact-match algorithms with pre-expanded variants
pub fn find_primer_expanded(
    cfg: &SearchConfig,
    read: &str,
    primer_variants: &[String],
) -> Option<PrimerMatch> {
    find_primer_bndm(cfg, read, primer_variants)
}
