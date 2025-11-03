use bio::alignment::pairwise::Aligner;
use anyhow;

/// Configuration for primer search
#[derive(Debug, Clone)]
pub struct SearchConfig {
    pub max_mismatches: usize,    // maximum allowed mismatches in primer sequence
    pub window: usize,            // number of bases to search from each end
    pub min_overlap: usize,       // minimum overlap length (like cutadapt -O)
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
        b'N' => true, // matches any base
        _ => false, // unknown code, treat as mismatch
    }
}

/// Find primer in read using semi-global alignment with IUPAC support
fn find_primer(cfg: &SearchConfig, read: &str, primer: &str) -> Option<PrimerMatch> {
    // Search in the first 'window' bases
    let window = cfg.window.min(read.len());
    let region = &read[..window];
    
    // Cutadapt-style scoring with IUPAC-aware matching
    let mut aligner = Aligner::with_capacity(
        primer.len(),
        region.len(),
        -1, // gap open
        -1, // gap extend
        |a, b| if matches_iupac(a, b) { 1 } else { -1 }, // IUPAC-aware
    );
    // Perform local alignment to find the best matching region (not just at position 0)
    let alignment = aligner.local(primer.as_bytes(), region.as_bytes());
    
    // Calculate coverage in both sequences
    let primer_coverage = alignment.xend - alignment.xstart;
    let read_coverage = alignment.yend - alignment.ystart;
    
    // Overlap is the maximum coverage between query (primer) and subject (read)
    let overlap_len = primer_coverage.max(read_coverage);
    
    // Print the overlap length for debugging
    println!("overlap_len: {}", overlap_len);
    
    // Check minimum overlap requirement
    if overlap_len < cfg.min_overlap.min(primer.len()) {
        return None;
    }
    
    // Count mismatches by manually checking with IUPAC awareness
    // Note: alignment.operations uses byte-level comparison, not our IUPAC-aware scoring
    // So we must manually iterate through aligned positions
    let compare_len = primer_coverage.min(read_coverage);
    let mut mismatches = 0;
    for i in 0..compare_len {
        if !matches_iupac(primer.as_bytes()[alignment.xstart + i], region.as_bytes()[alignment.ystart + i]) {
            mismatches += 1;
        }
    }
    
    // Debug: print the sequence being compared
    println!("  DEBUG: Aligned sequences - primer[{}..{}]='{}' vs region[{}..{}]='{}'", 
             alignment.xstart, alignment.xend,
             &primer[alignment.xstart..alignment.xend.min(primer.len())],
             alignment.ystart, alignment.yend,
             &region[alignment.ystart..alignment.yend.min(region.len())]);
    println!("  DEBUG: Alignment positions - primer starts at {}, region starts at {}, primer coverage={}, read coverage={}, alignment score={}", 
             alignment.xstart, alignment.ystart, primer_coverage, read_coverage, alignment.score);
    
    // Accept if mismatches are within allowed limit
    if mismatches <= cfg.max_mismatches {
        Some(PrimerMatch {
            start: alignment.ystart,
            end: alignment.yend,
            mismatches,
        })
    } else {
        println!("  → Rejected: {} mismatches (max allowed: {})", mismatches, cfg.max_mismatches);
        None
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
    
    println!("  DEBUG: Searching primer '{}' in end region (last 100bp): '{}'", 
             primer, if end_region.len() > 100 { &end_region[end_region.len()-100..] } else { end_region });
        
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
    println!("Trying Scenario 1: Forward primer at start...");
    if let Some(forward_match) = find_primer(cfg, read, forward_primer) {
        println!("  → Forward primer found at start (mismatches: {})", forward_match.mismatches);
        println!("  → Now searching for Rev-RC at end...");
        if let Some(reverse_match) = find_primer_at_end(cfg, read, &reverse_primer_rc, search_length) {
            println!("  → Rev-RC found at end (mismatches: {})", reverse_match.mismatches);
            println!(" ✅ Scenario 1 ACCEPTED: Forward orientation (Fwd at start, Rev-RC at end)");
            return PairedPrimerSearchResult {
                found: true,
                trim_start: forward_match.end,
                trim_end: reverse_match.start,
                needs_reverse_complement: false,
            };
        } else {
            println!("  ✗ Rev-RC NOT found at end");
        }
    } else {
        println!("  ✗ Forward primer NOT found at start");
    }
    
    // Scenario 2: Reverse primer at start, reverse complement of forward primer at end
    // Search in original read - trim positions are in original read coordinates
    // After trimming, the amplicon will be reverse complemented
    println!("Trying Scenario 2: Reverse primer at start...");
    if let Some(reverse_match) = find_primer(cfg, read, reverse_primer) {
        println!("  → Reverse primer found at start (mismatches: {})", reverse_match.mismatches);
        println!("  → Now searching for Fwd-RC at end...");
        if let Some(forward_match) = find_primer_at_end(cfg, read, &forward_primer_rc, search_length) {
            // Trim coordinates are in original read: trim from reverse_match.end to forward_match.start
            println!("  → Fwd-RC found at end (mismatches: {})", forward_match.mismatches);
            println!(" ✅ Scenario 2 ACCEPTED: Reverse orientation (Rev at start, Fwd-RC at end) - will reverse complement");
            return PairedPrimerSearchResult {
                found: true,
                trim_start: reverse_match.end,      // Start keeping after reverse primer
                trim_end: forward_match.start,      // Stop keeping before forward primer RC
                needs_reverse_complement: true,
            };
        } else {
            println!("  ✗ Fwd-RC NOT found at end");
        }
    } else {
        println!("  ✗ Reverse primer NOT found at start");
    }
    
    // Not found
    println!("❌ No primers found: Both scenarios failed");
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
        max_mismatches: usize,
        min_overlap: usize,
    ) -> anyhow::Result<Self> {
        Ok(Self {
            forward_primer,
            reverse_primer,
            search_length,
            config: SearchConfig {
                max_mismatches,
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
