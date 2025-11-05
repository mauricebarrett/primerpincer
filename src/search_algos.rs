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
//! Both algorithms (`find_primer_simd` and `find_primer_myers`) return the same type
//! and can be used interchangeably by `search_paired_primers` to locate primers at
//! the start and end of reads.

use crate::cli::Algorithm;
use crate::primer_search::SearchConfig;
use bio::alignment::distance::levenshtein;
use bio::alignment::{AlignmentOperation, pairwise::Aligner};
use bio::pattern_matching::myers::{Myers, MyersBuilder};
use once_cell::sync::Lazy;
use sassy::{Searcher, profiles::Iupac};

/// Result of a primer search
/// Contains the coordinates where a primer was found in a read sequence
#[derive(Debug, Clone)]
pub struct PrimerMatch {
    pub start: usize,
    pub end: usize,
}

/// Global MyersBuilder configured with IUPAC ambiguity codes
/// Created once and reused for all pattern building
static MYERS_BUILDER: Lazy<MyersBuilder> = Lazy::new(|| {
    let mut builder = MyersBuilder::new();

    // Configure IUPAC ambiguity codes
    let ambigs = vec![
        (b'M', b"AC".to_vec()),
        (b'R', b"AG".to_vec()),
        (b'W', b"AT".to_vec()),
        (b'S', b"CG".to_vec()),
        (b'Y', b"CT".to_vec()),
        (b'K', b"GT".to_vec()),
        (b'V', b"ACGMRS".to_vec()),
        (b'H', b"ACTMWY".to_vec()),
        (b'D', b"AGTRWK".to_vec()),
        (b'B', b"CGTSYK".to_vec()),
        (b'N', b"ACGTMRWSYKVHDB".to_vec()),
    ];

    for (base, equivalents) in ambigs {
        builder.ambig(base, equivalents.as_slice());
    }

    builder
});

/// Build a Myers pattern using the global MYERS_BUILDER
pub fn build_myers_pattern(primer: &str) -> Myers<u64> {
    MYERS_BUILDER.build_64(primer.as_bytes())
}

/// Find primer in read using a precompiled Myers pattern
pub fn find_primer_myers_precompiled(
    cfg: &SearchConfig,
    read: &str,
    myers: &mut Myers<u64>,
    primer_len: usize,
) -> Option<PrimerMatch> {
    // Search in the first 'window' bases
    let window = cfg.window.min(read.len());
    let region = &read[..window];
    let region_bytes = region.as_bytes();

    // Find all matches within the maximum edit distance
    let max_dist = cfg.edit_distance.min(u8::MAX as usize) as u8;
    let matches: Vec<_> = myers.find_all_lazy(region_bytes, max_dist).collect();

    // Find best match (lowest distance)
    let best_match =
        matches
            .iter()
            .min_by_key(|&&(_, dist)| dist)
            .and_then(|&(end_pos, _edit_distance)| {
                // Calculate start position
                let f_start = end_pos.saturating_sub(primer_len) + 1;
                let overlap_len = end_pos + 1 - f_start;

                // Check minimum overlap requirement
                let min_overlap_len = cfg.min_overlap.min(primer_len);
                if overlap_len >= min_overlap_len {
                    Some(PrimerMatch {
                        start: f_start,
                        end: end_pos,
                    })
                } else {
                    None
                }
            });

    best_match
}

/// Find primer in read using Myers algorithm (fallback for non-precompiled case)
fn find_primer_myers(cfg: &SearchConfig, read: &str, primer: &str) -> Option<PrimerMatch> {
    // Build pattern on the fly (slower, but works for backward compatibility)
    let mut myers = build_myers_pattern(primer);
    find_primer_myers_precompiled(cfg, read, &mut myers, primer.len())
}

/// Helper function to check if a base matches an IUPAC code
fn matches_iupac(iupac: u8, base: u8) -> bool {
    let iupac_upper = iupac.to_ascii_uppercase();
    let base_upper = base.to_ascii_uppercase();
    match iupac_upper {
        b'A' => base_upper == b'A',
        b'C' => base_upper == b'C',
        b'G' => base_upper == b'G',
        b'T' => base_upper == b'T',
        b'R' => base_upper == b'A' || base_upper == b'G', // A or G
        b'Y' => base_upper == b'C' || base_upper == b'T', // C or T
        b'M' => base_upper == b'A' || base_upper == b'C', // A or C
        b'K' => base_upper == b'G' || base_upper == b'T', // G or T
        b'S' => base_upper == b'C' || base_upper == b'G', // C or G
        b'W' => base_upper == b'A' || base_upper == b'T', // A or T
        b'B' => base_upper == b'C' || base_upper == b'G' || base_upper == b'T', // not A
        b'D' => base_upper == b'A' || base_upper == b'G' || base_upper == b'T', // not C
        b'H' => base_upper == b'A' || base_upper == b'C' || base_upper == b'T', // not G
        b'V' => base_upper == b'A' || base_upper == b'C' || base_upper == b'G', // not T
        b'N' => true,                                     // matches any base
        _ => false,                                       // unknown code, treat as mismatch
    }
}

fn find_primer_pairwise_local(cfg: &SearchConfig, read: &str, primer: &str) -> Option<PrimerMatch> {
    // Search in the first 'window' bases
    let window = cfg.window.min(read.len());
    let region = &read[..window];

    // Cutadapt-style scoring with IUPAC-aware matching
    let mut aligner = Aligner::with_capacity(
        primer.len(),
        region.len(),
        -1,                                              // gap open
        -1,                                              // gap extend
        |a, b| if matches_iupac(a, b) { 1 } else { -1 }, // IUPAC-aware
    );
    // Perform local alignment to find the best matching region (not just at position 0)
    let alignment = aligner.local(primer.as_bytes(), region.as_bytes());

    if alignment.operations.is_empty() {
        return None;
    }

    let mut mismatch_count = 0usize;
    for op in alignment.operations.iter() {
        match op {
            AlignmentOperation::Match => {}
            AlignmentOperation::Subst | AlignmentOperation::Ins | AlignmentOperation::Del => {
                mismatch_count += 1;
            }
            _ => {}
        }
    }

    if mismatch_count > cfg.max_mismatch {
        return None;
    }

    let start = alignment.ystart;
    let end = alignment.yend;
    if end <= start {
        return None;
    }

    let overlap_len = end - start;
    let min_required_overlap = cfg.min_overlap.min(primer.len());
    if overlap_len < min_required_overlap {
        return None;
    }

    Some(PrimerMatch { start, end })
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

    // Use sliding window approach with edit distance
    // Try all possible positions in the region
    for start in 0..=region.len().saturating_sub(min_overlap_len) {
        // Try different end positions to find the best match
        let max_end = (start + primer.len() + cfg.edit_distance).min(region.len());
        for end in (start + min_overlap_len)..=max_end {
            if end <= start {
                continue;
            }

            let region_slice = &region_bytes[start..end];

            // Calculate edit distance using standard Levenshtein
            let distance = levenshtein(primer_bytes, region_slice);
            let max_distance = cfg.edit_distance.min(u32::MAX as usize) as u32;

            if distance <= max_distance && distance < best_distance {
                let overlap_len = end - start;
                if overlap_len >= min_overlap_len {
                    best_match = Some(PrimerMatch { start, end });
                    best_distance = distance;
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

    let mut searcher = Searcher::<Iupac>::new_fwd();
    let matches = searcher.search(primer.as_bytes(), region_bytes, cfg.edit_distance);

    // Find best match (lowest cost)
    let best_match = matches.iter().min_by_key(|m| m.cost).and_then(|m| {
        let overlap_len = m.text_end - m.text_start;
        let min_overlap_len = cfg.min_overlap.min(primer.len());

        // Check minimum overlap requirement
        if overlap_len >= min_overlap_len {
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
        Algorithm::Myers => find_primer_myers(cfg, read, primer),
        Algorithm::Local => find_primer_pairwise_local(cfg, read, primer),
        Algorithm::Sassy => find_primer_sassy(cfg, read, primer),
    }
}
