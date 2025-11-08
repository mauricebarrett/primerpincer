use crate::cli::Algorithm;
use crate::preparing_input::PrimerSet;
use crate::search_algos::{PrimerMatch, find_primer, find_primer_bndm_cached, MyersPatternCache, BndmPatternCache};
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

/// Helper to search with optional cache using RefCell interior mutability
fn search_with_cache(
    cfg: &SearchConfig,
    read: &str,
    primer: &str,
    myers_cache: Option<&std::cell::RefCell<bio::pattern_matching::myers::long::Myers<u64>>>,
) -> Option<PrimerMatch> {
    if let Some(cache) = myers_cache {
        find_primer(cfg, read, primer, Some(&mut *cache.borrow_mut()))
    } else {
        find_primer(cfg, read, primer, None)
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
    myers_cache: Option<&MyersPatternCache>,
    bndm_cache: Option<&BndmPatternCache>,
) -> PairedPrimerSearchResult {
    // Helper to search with BNDM cache if available, otherwise use standard search
    let search_forward = |cfg: &SearchConfig, read: &str| -> Option<PrimerMatch> {
        if let Some(cache) = bndm_cache {
            find_primer_bndm_cached(cfg, read, &cache.forward)
        } else {
            search_with_cache(cfg, read, &primers.forward, myers_cache.map(|c| &c.forward))
        }
    };

    let search_reverse = |cfg: &SearchConfig, read: &str| -> Option<PrimerMatch> {
        if let Some(cache) = bndm_cache {
            find_primer_bndm_cached(cfg, read, &cache.reverse)
        } else {
            search_with_cache(cfg, read, &primers.reverse, myers_cache.map(|c| &c.reverse))
        }
    };

    let search_forward_rc = |cfg: &SearchConfig, read: &str| -> Option<PrimerMatch> {
        if let Some(cache) = bndm_cache {
            find_primer_bndm_cached(cfg, read, &cache.forward_rc)
        } else {
            search_with_cache(cfg, read, &primers.forward_rc, myers_cache.map(|c| &c.forward_rc))
        }
    };

    let search_reverse_rc = |cfg: &SearchConfig, read: &str| -> Option<PrimerMatch> {
        if let Some(cache) = bndm_cache {
            find_primer_bndm_cached(cfg, read, &cache.reverse_rc)
        } else {
            search_with_cache(cfg, read, &primers.reverse_rc, myers_cache.map(|c| &c.reverse_rc))
        }
    };

    // Scenario 1: Forward primer at start, reverse complement of reverse primer at end
    if let Some(forward_match) = search_forward(cfg, read) {
        if let Some(reverse_match) = {
            let search_len = search_length.min(read.len());
            let end_region = &read[read.len() - search_len..];
            if let Some(match_result) = search_reverse_rc(cfg, end_region) {
                let offset = read.len() - search_len;
                Some(PrimerMatch {
                    start: offset + match_result.start,
                    end: offset + match_result.end,
                })
            } else {
                None
            }
        } {
            return PairedPrimerSearchResult {
                found: true,
                trim_start: forward_match.end,
                trim_end: reverse_match.start,
                needs_reverse_complement: false,
            };
        }
    }

    // Scenario 2: Reverse primer at start, reverse complement of forward primer at end
    if let Some(reverse_match) = search_reverse(cfg, read) {
        if let Some(forward_match) = {
            let search_len = search_length.min(read.len());
            let end_region = &read[read.len() - search_len..];
            if let Some(match_result) = search_forward_rc(cfg, end_region) {
                let offset = read.len() - search_len;
                Some(PrimerMatch {
                    start: offset + match_result.start,
                    end: offset + match_result.end,
                })
            } else {
                None
            }
        } {
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

/// Matcher for finding and trimming primers from sequences
/// Cloned into each worker thread for private, lock-free pattern access
#[derive(Clone)]
pub struct PrimerMatcher {
    primers: PrimerSet,
    search_length: usize,
    config: SearchConfig,
    myers_cache: Option<MyersPatternCache>,
    bndm_cache: Option<BndmPatternCache>,
}

impl PrimerMatcher {
    /// Create a new PrimerMatcher with the given parameters
    pub fn new(
        primers: PrimerSet,
        search_length: usize,
        algorithm: Algorithm,
        error_rate: f64,
        min_overlap: usize,
        myers_cache: Option<MyersPatternCache>,
        bndm_cache: Option<BndmPatternCache>,
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
            myers_cache,
            bndm_cache,
        })
    }

    /// Search for paired primers in a sequence
    pub fn search_primers(&self, seq: &str) -> PairedPrimerSearchResult {
        search_paired_primers(&self.config, seq, &self.primers, self.search_length, self.myers_cache.as_ref(), self.bndm_cache.as_ref())
    }
}
