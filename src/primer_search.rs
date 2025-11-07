use crate::cli::Algorithm;
use crate::preparing_input::PrimerSet;
use crate::search_algos::{
    PrimerMatch, find_primer,
};
use anyhow;

/// Configuration for primer search
#[derive(Debug, Clone)]
pub struct SearchConfig {
    pub algorithm: Algorithm, // algorithm to use for matching
    pub error_rate: f64,      // maximum allowed error rate (errors / alignment_length)
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
}

impl PrimerMatcher {
    /// Create a new PrimerMatcher with the given parameters
    pub fn new(
        primers: PrimerSet,
        search_length: usize,
        algorithm: Algorithm,
        error_rate: f64,
        min_overlap: usize,
    ) -> anyhow::Result<Self> {
        Ok(Self {
            primers,
            search_length,
            config: SearchConfig {
                algorithm,
                error_rate,
                window: search_length,
                min_overlap,
            },
        })
    }

    /// Search for paired primers in a sequence
    pub fn search_primers(&self, seq: &str) -> PairedPrimerSearchResult {
        search_paired_primers(&self.config, seq, &self.primers, self.search_length)
    }
}
