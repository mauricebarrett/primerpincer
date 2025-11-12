//! Utilities for preparing primer input data

use once_cell::sync::Lazy;
use std::cell::RefCell;

use bio::pattern_matching::myers::{MyersBuilder, long::Myers};

/// Lookup table for reverse complement transformation
/// Maps each ASCII byte to its reverse complement
/// Works with both uppercase and lowercase input
static COMPLEMENT_TABLE: Lazy<[u8; 256]> = Lazy::new(|| {
    let mut table = [b'N'; 256];

    // Standard bases (uppercase)
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'N' as usize] = b'N';

    // IUPAC ambiguity codes (uppercase)
    table[b'R' as usize] = b'Y'; // A or G -> C or T
    table[b'Y' as usize] = b'R'; // C or T -> A or G
    table[b'M' as usize] = b'K'; // A or C -> G or T
    table[b'K' as usize] = b'M'; // G or T -> A or C
    table[b'S' as usize] = b'S'; // C or G -> C or G (symmetric)
    table[b'W' as usize] = b'W'; // A or T -> A or T (symmetric)
    table[b'B' as usize] = b'V'; // C or G or T (not A) -> A or C or G (not T)
    table[b'D' as usize] = b'H'; // A or G or T (not C) -> A or C or T (not G)
    table[b'H' as usize] = b'D'; // A or C or T (not G) -> A or G or T (not C)
    table[b'V' as usize] = b'B'; // A or C or G (not T) -> C or G or T (not A)

    // Lowercase versions
    table[b'a' as usize] = b'T';
    table[b't' as usize] = b'A';
    table[b'c' as usize] = b'G';
    table[b'g' as usize] = b'C';
    table[b'n' as usize] = b'N';
    table[b'r' as usize] = b'Y';
    table[b'y' as usize] = b'R';
    table[b'm' as usize] = b'K';
    table[b'k' as usize] = b'M';
    table[b's' as usize] = b'S';
    table[b'w' as usize] = b'W';
    table[b'b' as usize] = b'V';
    table[b'd' as usize] = b'H';
    table[b'h' as usize] = b'D';
    table[b'v' as usize] = b'B';

    table
});

/// Reverse complement a DNA sequence
/// Handles all IUPAC ambiguity codes correctly (uppercase and lowercase)
/// Optimized with lookup table: O(n) time, single pass in reverse
#[inline]
pub fn reverse_complement(seq: &str) -> String {
    let bytes = seq.as_bytes();
    let mut result = Vec::with_capacity(bytes.len());

    // Iterate backwards, using lookup table for instant mapping
    for &byte in bytes.iter().rev() {
        result.push(COMPLEMENT_TABLE[byte as usize]);
    }

    // Safe because we only produce valid ASCII from lookup table
    unsafe { String::from_utf8_unchecked(result) }
}

/// Prepared primer sequences with cached reverse complements
#[derive(Debug, Clone)]
pub struct PrimerSet {
    pub forward: String,
    pub reverse: String,
    pub forward_rc: String,
    pub reverse_rc: String,
}

impl PrimerSet {
    pub fn new(forward: impl Into<String>, reverse: impl Into<String>) -> Self {
        let forward = forward.into();
        let reverse = reverse.into();
        let forward_rc = reverse_complement(&forward);
        let reverse_rc = reverse_complement(&reverse);

        Self {
            forward,
            reverse,
            forward_rc,
            reverse_rc,
        }
    }
}

/// Expanded primer set with all concrete variants for exact-matching algorithms
/// Each primer may be degenerate, so we store Vec of concrete sequences
/// This is used by algorithms like BNDM that require concrete strings
#[derive(Debug, Clone)]
pub struct ExpandedPrimerSet {
    pub forward: Vec<String>,
    pub reverse: Vec<String>,
    pub forward_rc: Vec<String>,
    pub reverse_rc: Vec<String>,
}

impl ExpandedPrimerSet {
    /// Create a new expanded primer set by expanding degenerate primers into all concrete variants
    pub fn new(primers: &PrimerSet) -> Self {
        Self {
            forward: expand_degenerate_bases(&primers.forward),
            reverse: expand_degenerate_bases(&primers.reverse),
            forward_rc: expand_degenerate_bases(&primers.forward_rc),
            reverse_rc: expand_degenerate_bases(&primers.reverse_rc),
        }
    }
}

/// Myers pattern matchers prepared for all four primer orientations
#[derive(Clone)]
pub struct MyersPatternSet {
    pub forward: RefCell<Myers<u64>>,
    pub reverse: RefCell<Myers<u64>>,
    pub forward_rc: RefCell<Myers<u64>>,
    pub reverse_rc: RefCell<Myers<u64>>,
}

/// Global Myers builder with configured ambiguity tables
pub(crate) static MYERS_BUILDER: Lazy<MyersBuilder> = Lazy::new(|| {
    let mut builder = MyersBuilder::new();

    builder.ambig(b'M', b"AC");
    builder.ambig(b'R', b"AG");
    builder.ambig(b'W', b"AT");
    builder.ambig(b'S', b"CG");
    builder.ambig(b'Y', b"CT");
    builder.ambig(b'K', b"GT");
    builder.ambig(b'V', b"ACG");
    builder.ambig(b'H', b"ACT");
    builder.ambig(b'D', b"AGT");
    builder.ambig(b'B', b"CGT");
    builder.ambig(b'N', b"ACGT");

    builder.ambig(b'A', b"MRWVHDN");
    builder.ambig(b'C', b"MSYVHBN");
    builder.ambig(b'G', b"RSKDVBN");
    builder.ambig(b'T', b"WYKDHBN");

    builder
});

pub(crate) fn build_myers_matcher(pattern: &[u8]) -> Myers<u64> {
    MYERS_BUILDER.build_long(pattern)
}

impl MyersPatternSet {
    pub fn new(primers: &PrimerSet) -> Self {
        Self {
            forward: RefCell::new(build_myers_matcher(primers.forward.as_bytes())),
            reverse: RefCell::new(build_myers_matcher(primers.reverse.as_bytes())),
            forward_rc: RefCell::new(build_myers_matcher(primers.forward_rc.as_bytes())),
            reverse_rc: RefCell::new(build_myers_matcher(primers.reverse_rc.as_bytes())),
        }
    }
}

/// Expand IUPAC degenerative codes into all concrete base combinations
/// Returns all possible sequences represented by the degenerate code
///
/// # Example
/// "TCCTAGGGC" → ["TCCTAGGGC"] (no degeneracy)
/// "TCCTNGGGC" → ["TCCTAGGGC", "TCCTCGGGC", "TCCTGGGGC", "TCCTTGGGC"]
pub fn expand_degenerate_bases(seq: &str) -> Vec<String> {
    // IUPAC ambiguity code expansions
    let expand_base = |base: char| -> Vec<char> {
        match base.to_ascii_uppercase() {
            'A' => vec!['A'],
            'C' => vec!['C'],
            'G' => vec!['G'],
            'T' => vec!['T'],
            'R' => vec!['A', 'G'],           // puRine
            'Y' => vec!['C', 'T'],           // pYrimidine
            'M' => vec!['A', 'C'],           // aMino
            'K' => vec!['G', 'T'],           // Keto
            'S' => vec!['C', 'G'],           // Strong
            'W' => vec!['A', 'T'],           // Weak
            'B' => vec!['C', 'G', 'T'],      // not A
            'D' => vec!['A', 'G', 'T'],      // not C
            'H' => vec!['A', 'C', 'T'],      // not G
            'V' => vec!['A', 'C', 'G'],      // not T
            'N' => vec!['A', 'C', 'G', 'T'], // aNy base
            _ => vec![base],                 // unknown or already concrete
        }
    };

    let bases: Vec<Vec<char>> = seq.chars().map(|c| expand_base(c)).collect();

    // Generate cartesian product of all base options
    fn cartesian_product(
        bases: &[Vec<char>],
        index: usize,
        current: String,
        result: &mut Vec<String>,
    ) {
        if index == bases.len() {
            result.push(current);
            return;
        }

        for &base in &bases[index] {
            cartesian_product(bases, index + 1, format!("{}{}", current, base), result);
        }
    }

    let mut result = Vec::new();
    if !bases.is_empty() {
        cartesian_product(&bases, 0, String::new(), &mut result);
    }
    result
}
