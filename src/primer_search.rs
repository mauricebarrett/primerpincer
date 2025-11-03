use bio::pattern_matching::myers::MyersBuilder;
use bio::alignment::distance::simd;
use anyhow;
use crate::cli::Algorithm;

/// Configuration for primer search
#[derive(Debug, Clone)]
pub struct SearchConfig {
    pub algorithm: Algorithm,      // algorithm to use for matching
    pub edit_distance: usize,    // maximum allowed edit distance in primer sequence
    pub window: usize,            // number of bases to search from each end
    pub min_overlap: usize,       // minimum overlap length
}

/// Result of a primer search
#[derive(Debug, Clone)]
pub struct PrimerMatch {
    pub start: usize,
    pub end: usize,
    pub mismatches: usize,
}

/// Reverse complement a DNA sequence
/// Handles all IUPAC ambiguity codes correctly
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            // Standard bases
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            // IUPAC ambiguity codes
            'R' => 'Y',  // A or G -> C or T
            'Y' => 'R',  // C or T -> A or G
            'M' => 'K',  // A or C -> G or T
            'K' => 'M',  // G or T -> A or C
            'S' => 'S',  // C or G -> C or G (symmetric)
            'W' => 'W',  // A or T -> A or T (symmetric)
            'B' => 'V',  // C or G or T (not A) -> A or C or G (not T)
            'D' => 'H',  // A or G or T (not C) -> A or C or T (not G)
            'H' => 'D',  // A or C or T (not G) -> A or G or T (not C)
            'V' => 'B',  // A or C or G (not T) -> C or G or T (not A)
            _ => 'N',    // Unknown character -> N
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

/// Create a MyersBuilder configured with IUPAC ambiguity codes
fn create_myers_builder() -> MyersBuilder {
    let mut builder = MyersBuilder::new();
    
    // Configure IUPAC ambiguity codes
    // Using Vec<u8> to handle variable length ambigs
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
}

/// Find primer in read using SIMD-accelerated edit distance
fn find_primer_simd(cfg: &SearchConfig, read: &str, primer: &str) -> Option<PrimerMatch> {
    // Search in the first 'window' bases
    let window = cfg.window.min(read.len());
    let region = &read[..window];
    
    let primer_bytes = primer.as_bytes();
    let region_bytes = region.as_bytes();
    let min_overlap_len = cfg.min_overlap.min(primer.len());
    
    let mut best_match: Option<PrimerMatch> = None;
    let mut best_distance = u32::MAX;
    
    // Use sliding window approach with SIMD edit distance
    // Try all possible positions in the region
    for start in 0..=region.len().saturating_sub(min_overlap_len) {
        // Try different end positions to find the best match
        let max_end = (start + primer.len() + cfg.edit_distance).min(region.len());
        for end in (start + min_overlap_len)..=max_end {
            if end <= start {
                continue;
            }
            
            let region_slice = &region_bytes[start..end];
            
            // Calculate edit distance using SIMD (returns u32)
            let distance = simd::levenshtein(primer_bytes, region_slice);
            let edit_distance_u32 = cfg.edit_distance.min(u32::MAX as usize) as u32;
            
            if distance <= edit_distance_u32 && distance < best_distance {
                let overlap_len = end - start;
                if overlap_len >= min_overlap_len {
                    best_match = Some(PrimerMatch {
                        start,
                        end,
                        mismatches: distance as usize,
                    });
                    best_distance = distance;
                }
            }
        }
    }
    
    best_match
}

/// Find primer in read using Myers algorithm with IUPAC support
fn find_primer_myers(cfg: &SearchConfig, read: &str, primer: &str) -> Option<PrimerMatch> {
    // Search in the first 'window' bases
    let window = cfg.window.min(read.len());
    let region = &read[..window];
    
    // Convert to bytes for Myers algorithm
    let primer_bytes = primer.as_bytes();
    let region_bytes = region.as_bytes();
    
    // Create Myers matcher with IUPAC ambiguity support
    let builder = create_myers_builder();
    let mut myers = builder.build_64(primer_bytes);
    
    // Find all matches within the maximum edit distance
    // Convert to u8, clamping at u8::MAX
    let max_dist = cfg.edit_distance.min(u8::MAX as usize) as u8;
    let matches: Vec<_> = myers.find_all_lazy(region_bytes, max_dist).collect();
    
    // Find best match (lowest distance)
    let best_match = matches.iter()
        .min_by_key(|&&(_, dist)| dist)
        .and_then(|&(end_pos, edit_distance)| {
            // Calculate start position
            // Myers returns end positions (1-based in some contexts, but we use 0-based)
            // For approximate matches, start position is approximately end_pos - primer_len + 1
            let edit_dist = edit_distance as usize;
            println!("edit_dist: {}", edit_dist);
            let f_start = end_pos.saturating_sub(primer.len()) + 1;
            let overlap_len = end_pos + 1 - f_start;
            println!("overlap_len: {}", overlap_len);
            
            // Check minimum overlap requirement
            let min_overlap_len = cfg.min_overlap.min(primer.len());
            println!("min_overlap_len: {}", min_overlap_len);
            if overlap_len >= min_overlap_len {
                Some(PrimerMatch {
                    start: f_start,
                    end: end_pos,
                    mismatches: edit_dist,
                })
            } else {
                None
            }
        });
    
    best_match
}

/// Find primer in read using the selected algorithm
fn find_primer(cfg: &SearchConfig, read: &str, primer: &str) -> Option<PrimerMatch> {
    match cfg.algorithm {
        Algorithm::Simd => find_primer_simd(cfg, read, primer),
        Algorithm::Myers => find_primer_myers(cfg, read, primer),
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
    cfg: &SearchConfig,
    read: &str,
    forward_primer: &str,
    reverse_primer: &str,
    search_length: usize,
) -> PairedPrimerSearchResult {
    let forward_primer_rc = reverse_complement(forward_primer);
    let reverse_primer_rc = reverse_complement(reverse_primer);
    
    // Scenario 1: Forward primer at start, reverse complement of reverse primer at end
    println!("Trying Scenario 1: Forward primer at start, reverse complement of reverse primer at end");
    if let Some(forward_match) = find_primer(cfg, read, forward_primer) {
        if let Some(reverse_match) = find_primer_at_end(cfg, read, &reverse_primer_rc, search_length) {
            println!("✅ Scenario 1 PASSED: Found both primers");
            return PairedPrimerSearchResult {
                found: true,
                trim_start: forward_match.end,
                trim_end: reverse_match.start,
                needs_reverse_complement: false,
            };
        } else {
            println!("❌ Scenario 1 FAILED: Forward primer found, but reverse complement of reverse primer not found at end");
        }
    } else {
        println!("❌ Scenario 1 FAILED: Forward primer not found at start");
    }
    
    // Scenario 2: Reverse primer at start, reverse complement of forward primer at end
    // Search in original read - trim positions are in original read coordinates
    // After trimming, the amplicon will be reverse complemented
    println!("Trying Scenario 2: Reverse primer at start, reverse complement of forward primer at end");
    if let Some(reverse_match) = find_primer(cfg, read, reverse_primer) {
        if let Some(forward_match) = find_primer_at_end(cfg, read, &forward_primer_rc, search_length) {
            // Trim coordinates are in original read: trim from reverse_match.end to forward_match.start
            println!("✅ Scenario 2 PASSED: Found both primers (will reverse complement)");
            return PairedPrimerSearchResult {
                found: true,
                trim_start: reverse_match.end,      // Start keeping after reverse primer
                trim_end: forward_match.start,      // Stop keeping before forward primer RC
                needs_reverse_complement: true,
            };
        } else {
            println!("❌ Scenario 2 FAILED: Reverse primer found, but reverse complement of forward primer not found at end");
        }
    } else {
        println!("❌ Scenario 2 FAILED: Reverse primer not found at start");
    }
    
    // Not found
    println!("❌ Both scenarios FAILED: No primers found");
    PairedPrimerSearchResult {
        found: false,
        trim_start: 0,
        trim_end: 0,
        needs_reverse_complement: false,
    }
}

/// Matcher for finding and trimming primers from sequences
#[derive(Clone)]
pub struct PrimerMatcher {
    forward_primer: String,
    reverse_primer: String,
    search_length: usize,
    config: SearchConfig,
}

impl PrimerMatcher {
    /// Create a new PrimerMatcher with the given parameters
    pub fn new(
        forward_primer: String,
        reverse_primer: String,
        search_length: usize,
        algorithm: Algorithm,
        edit_distance: usize,
        min_overlap: usize,
    ) -> anyhow::Result<Self> {
        Ok(Self {
            forward_primer,
            reverse_primer,
            search_length,
            config: SearchConfig {
                algorithm,
                edit_distance,
                window: search_length,
                min_overlap,
            },
        })
    }
    
    /// Search for paired primers in a sequence
    pub fn search_primers(&self, seq: &str) -> PairedPrimerSearchResult {
        search_paired_primers(
            &self.config,
            seq,
            &self.forward_primer,
            &self.reverse_primer,
            self.search_length,
        )
    }
}
