use crate::preparing_input::PrimerSet;
use crate::preparing_input::{ExpandedPrimerSet, MyersPatternSet};
use crate::search_algos::Algorithm;
use crate::search_algos::{PrimerMatch, find_primer_degenerate, find_primer_expanded};
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

/// Search using degenerate-aware algorithms (Myers, Sassy)
fn search_with_degenerate(
    cfg: &SearchConfig,
    read: &str,
    primer: &str,
    myers_cache: Option<&std::cell::RefCell<bio::pattern_matching::myers::long::Myers<u64>>>,
) -> Option<PrimerMatch> {
    match cfg.algorithm {
        Algorithm::Myers => {
            if let Some(cache) = myers_cache {
                let mut borrowed = cache.borrow_mut();
                find_primer_degenerate(cfg, read, primer, Some(&mut *borrowed))
            } else {
                find_primer_degenerate(cfg, read, primer, None)
            }
        }
        Algorithm::Sassy => find_primer_degenerate(cfg, read, primer, None),
        Algorithm::Bndm => unreachable!("Exact-match algorithms should use search_with_expanded"),
        Algorithm::Hamming => {
            unreachable!("Exact-match algorithms should use search_with_expanded")
        }
    }
}

/// Search using exact-match algorithms with pre-expanded concrete variants
fn search_with_expanded(
    cfg: &SearchConfig,
    read: &str,
    expanded_variants: &[String],
) -> Option<PrimerMatch> {
    find_primer_expanded(cfg, read, expanded_variants)
}

/// Search for paired primers in a read
/// Scenario 1: Forward primer at start, reverse complement of reverse primer at end
/// Scenario 2: Reverse primer at start, reverse complement of forward primer at end (requires reverse complementing read)
pub fn search_paired_primers(
    cfg: &SearchConfig,
    read: &str,
    primers: &PrimerSet,
    search_length: usize,
    myers_patterns: Option<&MyersPatternSet>,
    expanded_primers: Option<&ExpandedPrimerSet>,
) -> PairedPrimerSearchResult {
    // Route to algorithm-specific search based on configuration
    match cfg.algorithm {
        Algorithm::Bndm | Algorithm::Hamming => search_paired_primers_expanded(
            cfg,
            read,
            primers,
            search_length,
            expanded_primers.expect("Exact-match algorithms require expanded primers"),
        ),
        _ => search_paired_primers_degenerate(cfg, read, primers, search_length, myers_patterns),
    }
}

/// Search for paired primers using Myers or Sassy (degenerate-aware)
fn search_paired_primers_degenerate(
    cfg: &SearchConfig,
    read: &str,
    primers: &PrimerSet,
    search_length: usize,
    myers_patterns: Option<&MyersPatternSet>,
) -> PairedPrimerSearchResult {
    // Scenario 1: Forward primer at start, reverse complement of reverse primer at end
    if let Some(forward_match) = search_with_degenerate(
        cfg,
        read,
        &primers.forward,
        myers_patterns.map(|c| &c.forward),
    ) {
        if let Some(reverse_match) = {
            let search_len = search_length.min(read.len());
            let end_region = &read[read.len() - search_len..];
            if let Some(match_result) = search_with_degenerate(
                cfg,
                end_region,
                &primers.reverse_rc,
                myers_patterns.map(|c| &c.reverse_rc),
            ) {
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
    if let Some(reverse_match) = search_with_degenerate(
        cfg,
        read,
        &primers.reverse,
        myers_patterns.map(|c| &c.reverse),
    ) {
        if let Some(forward_match) = {
            let search_len = search_length.min(read.len());
            let end_region = &read[read.len() - search_len..];
            if let Some(match_result) = search_with_degenerate(
                cfg,
                end_region,
                &primers.forward_rc,
                myers_patterns.map(|c| &c.forward_rc),
            ) {
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

/// Search for paired primers using exact-match algorithms with pre-expanded concrete sequences
fn search_paired_primers_expanded(
    cfg: &SearchConfig,
    read: &str,
    _primers: &PrimerSet,
    search_length: usize,
    expanded_primers: &ExpandedPrimerSet,
) -> PairedPrimerSearchResult {
    // Scenario 1: Forward primer at start, reverse complement of reverse primer at end
    if let Some(forward_match) = search_with_expanded(cfg, read, &expanded_primers.forward) {
        if let Some(reverse_match) = {
            let search_len = search_length.min(read.len());
            let end_region = &read[read.len() - search_len..];
            if let Some(match_result) =
                search_with_expanded(cfg, end_region, &expanded_primers.reverse_rc)
            {
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
    if let Some(reverse_match) = search_with_expanded(cfg, read, &expanded_primers.reverse) {
        if let Some(forward_match) = {
            let search_len = search_length.min(read.len());
            let end_region = &read[read.len() - search_len..];
            if let Some(match_result) =
                search_with_expanded(cfg, end_region, &expanded_primers.forward_rc)
            {
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
    myers_patterns: Option<MyersPatternSet>,
    expanded_primers: Option<ExpandedPrimerSet>,
}

impl PrimerMatcher {
    /// Create a new PrimerMatcher with the given parameters
    pub fn new(
        primers: PrimerSet,
        search_length: usize,
        algorithm: Algorithm,
        error_rate: f64,
        min_overlap: usize,
        myers_patterns: Option<MyersPatternSet>,
        expanded_primers: Option<ExpandedPrimerSet>,
    ) -> anyhow::Result<Self> {
        // Ensure expanded primers are available when BNDM or Hamming are selected
        let expanded_primers = match (algorithm, expanded_primers) {
            (Algorithm::Bndm | Algorithm::Hamming, Some(exp)) => Some(exp),
            (Algorithm::Bndm | Algorithm::Hamming, None) => Some(ExpandedPrimerSet::new(&primers)),
            (_, exp) => exp,
        };

        Ok(Self {
            primers,
            search_length,
            config: SearchConfig {
                algorithm,
                error_rate,
                window: search_length,
                min_overlap,
            },
            myers_patterns,
            expanded_primers,
        })
    }

    /// Search for paired primers in a sequence
    pub fn search_primers(&self, seq: &str) -> PairedPrimerSearchResult {
        search_paired_primers(
            &self.config,
            seq,
            &self.primers,
            self.search_length,
            self.myers_patterns.as_ref(),
            self.expanded_primers.as_ref(),
        )
    }
}
