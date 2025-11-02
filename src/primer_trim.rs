use crate::primer_search::{PrimerSearcher, HammingSearcher, AlignmentSearcher, SearchConfig, SearchMethod};

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

/// Enum to hold either searcher type
enum Searcher {
    Alignment(AlignmentSearcher),
    Hamming(HammingSearcher),
}

impl Searcher {
    fn find_primer(&self, read: &str, primer: &str) -> Option<crate::primer_search::PrimerMatch> {
        match self {
            Searcher::Alignment(s) => s.find_primer(read, primer),
            Searcher::Hamming(s) => s.find_primer(read, primer),
        }
    }
}

/// Matcher for finding and trimming primers from sequences
pub struct PrimerMatcher {
    forward_primer: String,
    forward_primer_rc: String,
    reverse_primer: String,
    reverse_primer_rc: String,
    searcher: Searcher,
}

impl PrimerMatcher {
    /// Create a new PrimerMatcher with the specified configuration
    pub fn new(
        forward_primer: String,
        reverse_primer: String,
        search_length: usize,
        max_mismatches: usize,
        algorithm: SearchMethod,
    ) -> Self {
        let searcher = match algorithm {
            SearchMethod::Alignment => {
                // For Alignment, max_error_rate defaults to 10% (0.1)
                let config = SearchConfig {
                    max_error_rate: 0.1,
                    max_mismatches, // Not used for Alignment but keep for consistency
                    window: search_length,
                    method: SearchMethod::Alignment,
                };
                Searcher::Alignment(AlignmentSearcher { cfg: config })
            }
            SearchMethod::Hamming => {
                let config = SearchConfig {
                    max_error_rate: 0.0, // Not used for Hamming
                    max_mismatches,
                    window: search_length,
                    method: SearchMethod::Hamming,
                };
                Searcher::Hamming(HammingSearcher { cfg: config })
            }
        };
        
        Self {
            forward_primer: forward_primer.clone(),
            forward_primer_rc: reverse_complement(&forward_primer),
            reverse_primer: reverse_primer.clone(),
            reverse_primer_rc: reverse_complement(&reverse_primer),
            searcher,
        }
    }
    
    /// Find any primer match in the first search_length bases
    /// Returns (primer_end_position, which_primer, orientation)
    /// - primer_end_position: position where the primer ends (start of trimmed sequence)
    /// - which_primer: true = forward primer, false = reverse primer
    /// - orientation: true = forward orientation, false = reverse complement
    pub fn find_primer_match(&self, seq: &str) -> Option<(usize, bool, bool)> {
        // Try forward primer forward orientation
        if let Some(m) = self.searcher.find_primer(seq, &self.forward_primer) {
            return Some((m.end, true, true));
        }
        // Try forward primer reverse complement
        if let Some(m) = self.searcher.find_primer(seq, &self.forward_primer_rc) {
            return Some((m.end, true, false));
        }
        // Try reverse primer forward orientation
        if let Some(m) = self.searcher.find_primer(seq, &self.reverse_primer) {
            return Some((m.end, false, true));
        }
        // Try reverse primer reverse complement
        if let Some(m) = self.searcher.find_primer(seq, &self.reverse_primer_rc) {
            return Some((m.end, false, false));
        }
        None
    }
    
    /// Trim sequence and quality based on primer match position
    /// Returns (trimmed_sequence, trimmed_quality) if trimming is successful, None otherwise
    pub fn trim_sequence(&self, seq: &str, qual: &str, trim_pos: usize) -> Option<(String, String)> {
        let trimmed_seq = &seq[trim_pos..];
        let trimmed_qual = if qual.len() > trim_pos {
            &qual[trim_pos..]
        } else {
            qual
        };
        
        // Don't keep sequences that become empty after trimming
        if trimmed_seq.is_empty() {
            return None;
        }
        
        Some((trimmed_seq.to_string(), trimmed_qual.to_string()))
    }
}
