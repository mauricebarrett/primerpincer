use bio::alignment::pairwise::*;
use anyhow;

/// Configuration for primer search
#[derive(Debug, Clone)]
pub struct SearchConfig {
    pub max_mismatches: usize,    // maximum allowed mismatches in primer sequence
    pub window: usize,            // number of bases to search from each end
}

/// Result of a primer search
#[derive(Debug, Clone)]
pub struct PrimerMatch {
    pub start: usize,
    pub end: usize,
    pub mismatches: usize,
}

/// Reverse complement a DNA sequence
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'N' | 'n' => 'N',
            _ => 'N',
        })
        .collect()
}

/// Result of a paired primer search (forward at start, reverse at end)
#[derive(Debug, Clone)]
pub struct PairedPrimerSearchResult {
    pub found: bool,
    pub trim_start: usize,  // Position to trim from start
    pub trim_end: usize,     // Position to trim from end (from 3' end)
    pub needs_reverse_complement: bool,  // Whether read needs to be reverse complemented
}

// Trait for primer search algorithms
pub(crate) trait PrimerSearcher {
    fn find_primer(&self, read: &str, primer: &str) -> Option<PrimerMatch>;
}

/// Search for primer at the end of a read (last search_length bases)
/// Searches in the last search_length bases for the primer
/// Returns the match with coordinates relative to the original read
fn find_primer_at_end(
    searcher: &dyn PrimerSearcher,
    read: &str,
    primer: &str,
    search_length: usize,
) -> Option<PrimerMatch> {
    if read.len() < search_length {
        return None;
    }
    let end_region = &read[read.len() - search_length..];
    
    // Search for primer in the end region
    if let Some(match_result) = searcher.find_primer(end_region, primer) {
        // Convert coordinates from end_region to original read
        let offset = read.len() - search_length;
        Some(PrimerMatch {
            start: offset + match_result.start,
            end: offset + match_result.end,
            mismatches: match_result.mismatches,
        })
    } else {
        None
    }
}

/// Search for paired primers in a read
/// Scenario 1: Forward primer at start, reverse complement of reverse primer at end
/// Scenario 2: Reverse primer at start, reverse complement of forward primer at end (requires reverse complementing read)
pub fn search_paired_primers(
    searcher: &dyn PrimerSearcher,
    read: &str,
    forward_primer: &str,
    reverse_primer: &str,
    search_length: usize,
) -> PairedPrimerSearchResult {
    let forward_primer_rc = reverse_complement(forward_primer);
    let reverse_primer_rc = reverse_complement(reverse_primer);
    
    // Scenario 1: Forward primer at start, reverse complement of reverse primer at end
    if let Some(forward_match) = searcher.find_primer(read, forward_primer) {
        if let Some(reverse_match) = find_primer_at_end(searcher, read, &reverse_primer_rc, search_length) {
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
    if let Some(reverse_match) = searcher.find_primer(read, reverse_primer) {
        if let Some(forward_match) = find_primer_at_end(searcher, read, &forward_primer_rc, search_length) {
            // Trim coordinates are in original read: trim from reverse_match.end to forward_match.start
            return PairedPrimerSearchResult {
                found: true,
                trim_start: reverse_match.end,      // Start keeping after reverse primer
                trim_end: forward_match.start,      // Stop keeping before forward primer RC
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

//
// ---------------- ALIGNMENT (semi-global alignment) ----------------
//
#[derive(Clone)]
pub struct AlignmentSearcher {
    pub cfg: SearchConfig,
}

impl PrimerSearcher for AlignmentSearcher {
    fn find_primer(&self, read: &str, primer: &str) -> Option<PrimerMatch> {
        let window = self.cfg.window.min(read.len());
        let region = &read[..window];
        let mut aligner = Aligner::with_capacity(
            primer.len(),
            region.len(),
            -5, // gap open
            -1, // gap extend
            |a, b| if a == b { 2 } else { -3 },
        );
        let alignment = aligner.semiglobal(primer.as_bytes(), region.as_bytes());

        // Compute edit distance from alignment result
        let aligned_primer_len = alignment.xend - alignment.xstart;
        let aligned_region_len = alignment.yend - alignment.ystart;
        
        let mut errors = 0;
        let primer_bytes = primer.as_bytes();
        let region_bytes = region.as_bytes();
        
        // Count mismatches in the aligned portion
        let min_len = aligned_primer_len.min(aligned_region_len);
        for i in 0..min_len {
            let primer_idx = alignment.xstart + i;
            let region_idx = alignment.ystart + i;
            if primer_idx < primer_bytes.len() && region_idx < region_bytes.len() {
                if primer_bytes[primer_idx] != region_bytes[region_idx] {
                    errors += 1;
                }
            }
        }
        
        // Add penalty for length differences (indels)
        errors += aligned_primer_len.abs_diff(aligned_region_len);
        
        if errors <= self.cfg.max_mismatches {
            Some(PrimerMatch {
                start: alignment.ystart,
                end: alignment.yend,
                mismatches: errors,
            })
        } else {
            None
        }
    }
}


/// Matcher for finding and trimming primers from sequences
#[derive(Clone)]
pub struct PrimerMatcher {
    forward_primer: String,
    reverse_primer: String,
    search_length: usize,
    searcher: AlignmentSearcher,
}

impl PrimerMatcher {
    /// Create a new PrimerMatcher with the specified configuration
    /// Uses semi-global alignment for primer matching
    pub fn new(
        forward_primer: String,
        reverse_primer: String,
        search_length: usize,
        max_mismatches: usize,
    ) -> anyhow::Result<Self> {
        // Use semi-global alignment with max mismatches
        let config = SearchConfig {
            max_mismatches,
            window: search_length,
        };
        let searcher = AlignmentSearcher { cfg: config };
        
        Ok(Self {
            forward_primer,
            reverse_primer,
            search_length,
            searcher,
        })
    }
    
    /// Search for paired primers and return result
    pub fn search_primers(&self, seq: &str) -> PairedPrimerSearchResult {
        search_paired_primers(
            &self.searcher,
            seq,
            &self.forward_primer,
            &self.reverse_primer,
            self.search_length,
        )
    }
}
