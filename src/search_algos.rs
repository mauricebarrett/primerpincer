//! Search algorithms for finding primers in DNA sequences.
//!
//! This module contains implementations of different search algorithms for finding
//! primer sequences within reads. Each search algorithm returns an `Option<PrimerMatch>`
//! to `search_paired_primers` where:
//! - `Some(PrimerMatch)` indicates a match was found with the best match details.
//!   The `PrimerMatch` contains `start` and `end` coordinates (0-based) that are used
//!   by `search_paired_primers` to determine trim positions for removing primers from reads.
//!   The coordinates indicate where the primer was found in the read sequence.
//! - `None` indicates no match was found within the specified constraints.
//!
//! The search algorithms return the same type (`PrimerMatch`) and can be used
//! by `search_paired_primers` to locate primers at the start and end of reads.

use crate::cli::Algorithm;
use crate::primer_search::SearchConfig;
use bio::alignment::distance::levenshtein;
use sassy::{Searcher, profiles::Iupac};

/// Result of a primer search
/// Contains the coordinates where a primer was found in a read sequence
#[derive(Debug, Clone)]
pub struct PrimerMatch {
    pub start: usize,
    pub end: usize,
}

/// Find primer in read using Levenshtein edit distance
fn find_primer_levenshtein(cfg: &SearchConfig, read: &str, primer: &str) -> Option<PrimerMatch> {
    // Search in the first 'window' bases
    let window = cfg.window.min(read.len());
    let region = &read[..window];

    let primer_bytes = primer.as_bytes();
    let region_bytes = region.as_bytes();
    let min_overlap_len = cfg.min_overlap.min(primer.len());

    let mut best_match: Option<PrimerMatch> = None;
    let mut best_distance = u32::MAX;

    // Use sliding window approach with error rate threshold
    // Try all possible positions in the region
    for start in 0..=region.len().saturating_sub(min_overlap_len) {
        // Try different end positions to find the best match
        // Use generous upper bound for end positions
        let max_end = (start + (primer.len() as f64 / (1.0 - cfg.error_rate)).ceil() as usize).min(region.len());
        for end in (start + min_overlap_len)..=max_end {
            if end <= start {
                continue;
            }

            let region_slice = &region_bytes[start..end];
            let overlap_len = end - start;

            // Calculate edit distance using standard Levenshtein
            let distance = levenshtein(primer_bytes, region_slice);
            let error_rate = (distance as f64) / (overlap_len as f64);

            if error_rate <= cfg.error_rate && (distance as u32) < best_distance {
                if overlap_len >= min_overlap_len {
                    best_match = Some(PrimerMatch { start, end });
                    best_distance = distance as u32;
                }
            }
        }
    }

    best_match
}

/// Find primer in read using Sassy SIMD-accelerated search
fn find_primer_sassy(cfg: &SearchConfig, read: &str, primer: &str) -> Option<PrimerMatch> {
    // Search in the first 'window' bases
    let window = cfg.window.min(read.len());
    let region_bytes = &read[..window].as_bytes();

    // Use overhang cost to handle partial primer matches at read boundaries
    // Overhang cost of 0.5 means each overhanging character is penalized by 0.5
    // instead of the full cost - important for truncated primers at contig/read ends
    let overhang_cost = 0.0;
    let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(overhang_cost);
    // Convert error_rate to max allowed edits (generous upper bound)
    let max_edits = ((primer.len() as f64) * cfg.error_rate / (1.0 - cfg.error_rate)).ceil() as usize;
    let matches = searcher.search(primer.as_bytes(), region_bytes, max_edits);

    // Find best match (lowest cost)
    let best_match = matches.iter().min_by_key(|m| m.cost).and_then(|m| {
        let overlap_len = m.text_end - m.text_start;
        let min_overlap_len = cfg.min_overlap.min(primer.len());
        let error_rate = (m.cost as f64) / (overlap_len as f64);
        
        // Extract the aligned region from the read
        let aligned_region = std::str::from_utf8(&region_bytes[m.text_start..m.text_end])
            .unwrap_or("<invalid utf8>");
        
        // Debug: print match details with aligned region and error rate
        eprintln!("  [Sassy] Match: pos=[{}..{}] len={} primer_len={} cost={} error_rate={:.2}%", 
                 m.text_start, m.text_end, overlap_len, primer.len(), m.cost, error_rate * 100.0);
        eprintln!("    Primer:  {}", std::str::from_utf8(primer.as_bytes()).unwrap_or("<invalid>"));
        eprintln!("    Aligned: {}", aligned_region);

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

/// Find primer in read using the selected algorithm
pub fn find_primer(cfg: &SearchConfig, read: &str, primer: &str) -> Option<PrimerMatch> {
    match cfg.algorithm {
        Algorithm::Levenshtein => find_primer_levenshtein(cfg, read, primer),
        Algorithm::Sassy => find_primer_sassy(cfg, read, primer),
    }
}
