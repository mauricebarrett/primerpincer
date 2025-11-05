use crate::cli::Algorithm;
use crate::preparing_input::PrimerSet;
use crate::search_algos::{
    PrimerMatch, build_myers_pattern, find_primer, find_primer_myers_precompiled,
};
use anyhow;
use bio::pattern_matching::myers::Myers;

/// Configuration for primer search
#[derive(Debug, Clone)]
pub struct SearchConfig {
    pub algorithm: Algorithm, // algorithm to use for matching
    pub edit_distance: usize, // maximum allowed edit distance in primer sequence
    pub max_mismatch: usize,  // maximum allowed mismatches for local alignment
    pub window: usize,        // number of bases to search from each end
    pub min_overlap: usize,   // minimum overlap length
}

/// Result of a paired primer search (forward at start, reverse at end)
#[derive(Debug, Clone)]
pub struct PairedPrimerSearchResult {
    pub found: bool,
    pub trim_start: usize,              // Position to trim from start
    pub trim_end: usize,                // Position to trim from end (from 3' end)
    pub needs_reverse_complement: bool, // Whether read needs to be reverse complemented
}

/// Precompiled Myers patterns for all 4 primers
/// Designed to be cloned into each worker thread for zero-lock access
#[derive(Clone)]
pub struct PrecompiledMyersPatterns {
    pub forward: Myers<u64>,
    pub reverse: Myers<u64>,
    pub forward_rc: Myers<u64>,
    pub reverse_rc: Myers<u64>,
    pub forward_len: usize,
    pub reverse_len: usize,
    pub forward_rc_len: usize,
    pub reverse_rc_len: usize,
}

impl PrecompiledMyersPatterns {
    /// Build precompiled Myers patterns from a PrimerSet
    pub fn from_primer_set(primers: &PrimerSet) -> Self {
        Self {
            forward: build_myers_pattern(&primers.forward),
            reverse: build_myers_pattern(&primers.reverse),
            forward_rc: build_myers_pattern(&primers.forward_rc),
            reverse_rc: build_myers_pattern(&primers.reverse_rc),
            forward_len: primers.forward.len(),
            reverse_len: primers.reverse.len(),
            forward_rc_len: primers.forward_rc.len(),
            reverse_rc_len: primers.reverse_rc.len(),
        }
    }
}

/// Search for primer at the end of a read (last search_length bases)
/// Returns the match with coordinates relative to the original read
fn find_primer_at_end(
    cfg: &SearchConfig,
    read: &str,
    primer: &str,
    search_length: usize,
) -> Option<PrimerMatch> {
    // Find the start position of the search region
    let search_len = search_length.min(read.len());
    // Extract the end region (last search_length bases)
    let end_region = &read[read.len() - search_len..];

    // Perform semi-global alignment with IUPAC support
    if let Some(match_result) = find_primer(cfg, end_region, primer) {
        // Convert coordinates from end_region to original read
        let offset = read.len() - search_len;
        Some(PrimerMatch {
            start: offset + match_result.start,
            end: offset + match_result.end,
        })
    } else {
        None
    }
}

/// Search for primer at the end using precompiled Myers pattern
fn find_primer_at_end_myers(
    cfg: &SearchConfig,
    read: &str,
    myers: &mut Myers<u64>,
    primer_len: usize,
    search_length: usize,
) -> Option<PrimerMatch> {
    let search_len = search_length.min(read.len());
    let end_region = &read[read.len() - search_len..];

    if let Some(match_result) = find_primer_myers_precompiled(cfg, end_region, myers, primer_len) {
        let offset = read.len() - search_len;
        Some(PrimerMatch {
            start: offset + match_result.start,
            end: offset + match_result.end,
        })
    } else {
        None
    }
}

/// Search for paired primers using precompiled Myers patterns
/// Each thread has its own mutable copy of patterns (no cloning needed)
fn search_paired_primers_myers(
    cfg: &SearchConfig,
    read: &str,
    patterns: &mut PrecompiledMyersPatterns,
    search_length: usize,
) -> PairedPrimerSearchResult {
    // Scenario 1: Forward primer at start, reverse complement of reverse primer at end
    if let Some(forward_match) =
        find_primer_myers_precompiled(cfg, read, &mut patterns.forward, patterns.forward_len)
    {
        if let Some(reverse_match) = find_primer_at_end_myers(
            cfg,
            read,
            &mut patterns.reverse_rc,
            patterns.reverse_rc_len,
            search_length,
        ) {
            return PairedPrimerSearchResult {
                found: true,
                trim_start: forward_match.end,
                trim_end: reverse_match.start,
                needs_reverse_complement: false,
            };
        }
    }

    // Scenario 2: Reverse primer at start, reverse complement of forward primer at end
    if let Some(reverse_match) =
        find_primer_myers_precompiled(cfg, read, &mut patterns.reverse, patterns.reverse_len)
    {
        if let Some(forward_match) = find_primer_at_end_myers(
            cfg,
            read,
            &mut patterns.forward_rc,
            patterns.forward_rc_len,
            search_length,
        ) {
            return PairedPrimerSearchResult {
                found: true,
                trim_start: reverse_match.end,
                trim_end: forward_match.start,
                needs_reverse_complement: true,
            };
        }
    }

    // Not found
    PairedPrimerSearchResult {
        found: false,
        trim_start: 0,
        trim_end: 0,
        needs_reverse_complement: false,
    }
}

/// Search for paired primers in a read
/// Scenario 1: Forward primer at start, reverse complement of reverse primer at end
/// Scenario 2: Reverse primer at start, reverse complement of forward primer at end (requires reverse complementing read)
pub fn search_paired_primers(
    cfg: &SearchConfig,
    read: &str,
    primers: &PrimerSet,
    search_length: usize,
) -> PairedPrimerSearchResult {
    // Scenario 1: Forward primer at start, reverse complement of reverse primer at end
    if let Some(forward_match) = find_primer(cfg, read, &primers.forward) {
        if let Some(reverse_match) =
            find_primer_at_end(cfg, read, &primers.reverse_rc, search_length)
        {
            return PairedPrimerSearchResult {
                found: true,
                trim_start: forward_match.end,
                trim_end: reverse_match.start,
                needs_reverse_complement: false,
            };
        }
    }

    // Scenario 2: Reverse primer at start, reverse complement of forward primer at end
    // Search in original read - trim positions are in original read coordinates
    // After trimming, the amplicon will be reverse complemented
    if let Some(reverse_match) = find_primer(cfg, read, &primers.reverse) {
        if let Some(forward_match) =
            find_primer_at_end(cfg, read, &primers.forward_rc, search_length)
        {
            // Trim coordinates are in original read: trim from reverse_match.end to forward_match.start
            return PairedPrimerSearchResult {
                found: true,
                trim_start: reverse_match.end, // Start keeping after reverse primer
                trim_end: forward_match.start, // Stop keeping before forward primer RC
                needs_reverse_complement: true,
            };
        }
    }

    // Not found
    PairedPrimerSearchResult {
        found: false,
        trim_start: 0,
        trim_end: 0,
        needs_reverse_complement: false,
    }
}

/// Matcher for finding and trimming primers from sequences
/// Cloned into each worker thread for private, lock-free pattern access
#[derive(Clone)]
pub struct PrimerMatcher {
    primers: PrimerSet,
    search_length: usize,
    config: SearchConfig,
    /// Patterns are cloned into each thread - no locks needed
    myers_patterns: Option<PrecompiledMyersPatterns>,
}

impl PrimerMatcher {
    /// Create a new PrimerMatcher with the given parameters
    pub fn new(
        primers: PrimerSet,
        search_length: usize,
        algorithm: Algorithm,
        edit_distance: usize,
        max_mismatch: usize,
        min_overlap: usize,
        myers_patterns: Option<PrecompiledMyersPatterns>,
    ) -> anyhow::Result<Self> {
        // Use provided patterns or build them on the fly (for backward compatibility)
        let myers_patterns = myers_patterns.or_else(|| {
            if matches!(algorithm, Algorithm::Myers) {
                Some(PrecompiledMyersPatterns::from_primer_set(&primers))
            } else {
                None
            }
        });

        Ok(Self {
            primers,
            search_length,
            config: SearchConfig {
                algorithm,
                edit_distance,
                max_mismatch,
                window: search_length,
                min_overlap,
            },
            myers_patterns,
        })
    }

    /// Search for paired primers in a sequence
    /// Takes &mut self to allow mutable access to Myers patterns (no cloning needed)
    pub fn search_primers(&mut self, seq: &str) -> PairedPrimerSearchResult {
        // Use precompiled Myers patterns if available
        if let Some(ref mut patterns) = self.myers_patterns {
            search_paired_primers_myers(&self.config, seq, patterns, self.search_length)
        } else {
            search_paired_primers(&self.config, seq, &self.primers, self.search_length)
        }
    }
}
