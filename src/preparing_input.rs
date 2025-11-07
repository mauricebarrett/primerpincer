//! Utilities for preparing primer input data

use once_cell::sync::Lazy;

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
