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
#[derive(ValueEnum, Clone, Debug, Copy, Default)]
pub enum Algorithm {
    /// Pattern matching algorithm as described in Beeloo and Koerkamp (2025)
    #[default]
    Sassy,
    /// Rust Bio's Myers bit-parallel approximate pattern matching algorithm as described in Myers (1999). Implementation is very similar to Edlib's (Šošić and Šikić, 2017).
    Myers,
    /// Hamming distance algorithm as described in Waterman and Eggert (1987). Can tolerate mismatches but not indels.
    Hamming,
    /// Rust Bio's BNDM exact pattern matching algorithm as described in Baeza-Yates and Gonnet (1992). Exact matching only. No mismatch or indels tolerated.
    Bndm,
}
use crate::preparing_input::build_myers_matcher;
use crate::primer_search::SearchConfig;
use bio::alignment::distance::simd::hamming as simd_hamming;
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

/// Find primer in read using Sassy SIMD-accelerated search
fn find_primer_sassy(cfg: &SearchConfig, read: &str, primer: &str) -> Option<PrimerMatch> {
    // Search in the first 'window' bases
    let window = cfg.window.min(read.len());
    let region_bytes = &read.as_bytes()[..window];

    // Use overhang cost to handle partial primer matches at read boundaries
    // Overhang cost of 0.5 means each overhanging character is penalized by 0.5
    // instead of the full cost - important for truncated primers at contig/read ends
    let overhang_cost = 0.5;
    let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(overhang_cost);
    // Convert error_rate to max allowed edits (generous upper bound)
    let max_edits =
        ((primer.len() as f64) * cfg.error_rate / (1.0 - cfg.error_rate)).ceil() as usize;
    let matches = searcher.search(primer.as_bytes(), &region_bytes, max_edits);

    // Find best match (lowest cost)
    matches.iter().min_by_key(|m| m.cost).and_then(|m| {
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
    })
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
    let region_bytes = &read.as_bytes()[..window];
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

/// Find primer in read using Hamming distance algorithm
fn find_primer_hamming(
    cfg: &SearchConfig,
    read: &str,
    primer_variants: &[String],
) -> Option<PrimerMatch> {
    let window = cfg.window.min(read.len());
    let region_bytes = &read.as_bytes()[..window];

    // Guard: no variants provided
    if primer_variants.is_empty() {
        return None;
    }

    let min_overlap_len = cfg.min_overlap.min(primer_variants[0].len());

    // Track best match across all variants and positions
    let mut best_match: Option<PrimerMatch> = None;
    let mut best_distance: u64 = u64::MAX;

    for primer_variant in primer_variants {
        let primer_bytes = primer_variant.as_bytes();
        let primer_len = primer_bytes.len();

        // Slide the primer across the searchable region (no indels; compare fixed-length windows)
        for start in 0..window {
            let max_overlap = window - start;
            if max_overlap == 0 {
                break;
            }
            let overlap_len = primer_len.min(max_overlap);
            if overlap_len < min_overlap_len {
                continue;
            }

            // Compute Hamming distance on the overlapping portion
            let p_slice = &primer_bytes[..overlap_len];
            let r_slice = &region_bytes[start..start + overlap_len];
            let distance: u64 = simd_hamming(p_slice, r_slice);

            // Check error rate threshold
            let error_rate = (distance as f64) / (overlap_len as f64);
            if error_rate > cfg.error_rate {
                continue;
            }

            if distance < best_distance {
                best_distance = distance;
                best_match = Some(PrimerMatch {
                    start,
                    end: start + overlap_len - 1, // inclusive end
                });

                // Early exit on perfect match
                if best_distance == 0 {
                    return best_match;
                }
            }
        }
    }

    best_match
}

// Find primer in read using BNDM for exact matching
/// Searches all concrete variants (expanded from degenerate codes)
/// Note: BNDM requires complete exact matches of the full pattern, so min_overlap
/// is not applicable (overlap_len always equals primer length for any match)
fn find_primer_bndm(
    cfg: &SearchConfig,
    read: &str,
    primer_variants: &[String],
) -> Option<PrimerMatch> {
    // Search in the first 'window' bases
    let window = cfg.window.min(read.len());
    let region_bytes = &read.as_bytes()[..window];

    // Try each primer variant (from degenerate expansion)
    for primer_variant in primer_variants {
        let primer_bytes = primer_variant.as_bytes();
        let matcher = BNDM::new(primer_bytes);

        // Find first exact match in the window
        // BNDM always matches the complete pattern, so no partial overlap check needed
        if let Some(pos) = matcher.find_all(region_bytes).next() {
            let match_end = pos + primer_bytes.len();
            return Some(PrimerMatch {
                start: pos,
                end: match_end - 1, // Make end inclusive to match Myers convention
            });
        }
    }

    None
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
        Algorithm::Bndm => unreachable!("Requires pre-expanded variants"),
        Algorithm::Hamming => unreachable!("Requires pre-expanded variants"),
    }
}

/// Find primer in read using exact-match algorithms with pre-expanded variants
pub fn find_primer_expanded(
    cfg: &SearchConfig,
    read: &str,
    primer_variants: &[String],
) -> Option<PrimerMatch> {
    match cfg.algorithm {
        Algorithm::Bndm => find_primer_bndm(cfg, read, primer_variants),
        Algorithm::Hamming => find_primer_hamming(cfg, read, primer_variants),
        _ => unreachable!("Degenerate-aware algorithms should use find_primer_degenerate"),
    }
}
