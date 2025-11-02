use bio::alignment::pairwise::*;
use bio::alignment::AlignmentMode;

/// Configuration for primer search
#[derive(Debug, Clone)]
pub struct SearchConfig {
    pub max_error_rate: f32,     // e.g. 0.1 for 10%
    pub max_mismatches: usize,   // used by hamming mode
    pub window: usize,           // number of bases to search from each end
    pub method: SearchMethod,
}

/// Available algorithms
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SearchMethod {
    Alignment, // semi-global alignment
    Hamming,   // Hamming-distance sliding window
}

/// Result of a primer search
#[derive(Debug, Clone)]
pub struct PrimerMatch {
    pub start: usize,
    pub end: usize,
    pub mismatches: usize,
}

/// Trait that both algorithms implement
pub trait PrimerSearcher {
    fn find_primer(&self, read: &str, primer: &str) -> Option<PrimerMatch>;
}

//
// ---------------- ALIGNMENT (semi-global alignment) ----------------
//
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

        let errors = alignment.edit_distance;
        let max_errors = (primer.len() as f32 * self.cfg.max_error_rate).ceil() as usize;
        if errors <= max_errors {
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

//
// ---------------- HAMMING (Hamming distance) ----------------
//
pub struct HammingSearcher {
    pub cfg: SearchConfig,
}

impl PrimerSearcher for HammingSearcher {
    fn find_primer(&self, read: &str, primer: &str) -> Option<PrimerMatch> {
        let window = self.cfg.window.min(read.len());
        let region = &read[..window];
        let primer_len = primer.len();
        for i in 0..=region.len().saturating_sub(primer_len) {
            let mismatches = primer
                .bytes()
                .zip(region[i..].bytes())
                .filter(|(a, b)| a != b)
                .count();
            if mismatches <= self.cfg.max_mismatches {
                return Some(PrimerMatch {
                    start: i,
                    end: i + primer_len,
                    mismatches,
                });
            }
        }
        None
    }
}
