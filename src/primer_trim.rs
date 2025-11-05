use crate::cli::Algorithm;
use crate::preparing_input::{PrimerSet, reverse_complement};
use crate::primer_search::{PairedPrimerSearchResult, PrecompiledMyersPatterns, PrimerMatcher};
use anyhow;
use paraseq::Record;
use paraseq::parallel::{IntoProcessError, ParallelProcessor};
use std::borrow::Cow;
use std::io::Write;
use std::sync::Arc;
use std::sync::Mutex;

/// Trim sequence and quality based on search result
/// Handles reverse complementing if needed (after trimming)
/// Returns (trimmed_sequence, trimmed_quality) if trimming is successful, None otherwise
/// Uses Cow<str> to avoid allocations when reverse complementing isn't needed
fn trim_sequence<'a>(
    seq: &'a str,
    qual: &'a str,
    result: &PairedPrimerSearchResult,
) -> Option<(Cow<'a, str>, Cow<'a, str>)> {
    if !result.found {
        return None;
    }

    // Trim from both ends first (using original read coordinates)
    // trim_start: position to start keeping (trim everything before this)
    // trim_end: position to stop keeping (trim everything from this position to end)
    if result.trim_end <= result.trim_start {
        // Invalid trimming coordinates
        return None;
    }

    // Ensure we have valid trim coordinates
    let trim_end = result.trim_end.min(seq.len());
    if trim_end <= result.trim_start {
        return None; // Invalid or empty result
    }

    // Trim the sequence and quality
    let trimmed_seq = &seq[result.trim_start..trim_end];
    let trimmed_qual = qual.get(result.trim_start..trim_end)?; // Returns None if out of bounds

    // Reverse complement the trimmed amplicon if needed
    if result.needs_reverse_complement {
        let seq_rc = reverse_complement(trimmed_seq);
        // Optimize quality reversal: FASTQ quality is ASCII-only, so use byte operations
        // This is faster than .chars().rev() which requires UTF-8 decoding
        let mut qual_bytes = trimmed_qual.as_bytes().to_vec();
        qual_bytes.reverse();
        // Safe because FASTQ quality is always valid ASCII
        let qual_rc = unsafe { String::from_utf8_unchecked(qual_bytes) };
        Some((Cow::Owned(seq_rc), Cow::Owned(qual_rc)))
    } else {
        // Zero-copy: just borrow the slices without allocating
        Some((Cow::Borrowed(trimmed_seq), Cow::Borrowed(trimmed_qual)))
    }
}

/// Processor for trimming primers from FASTQ records
/// Implements ParallelProcessor for use with paraseq parallel processing
#[derive(Clone)]
pub struct PrimerTrimmer {
    matcher: PrimerMatcher,
    local_output: String,
    global_output: Arc<Mutex<Box<dyn Write + Send>>>,
}

impl PrimerTrimmer {
    pub fn new(
        output: Box<dyn Write + Send>,
        primers: PrimerSet,
        search_length: usize,
        algorithm: Algorithm,
        edit_distance: usize,
        max_mismatch: usize,
        min_overlap: usize,
        myers_patterns: Option<PrecompiledMyersPatterns>,
    ) -> anyhow::Result<Self> {
        let matcher = PrimerMatcher::new(
            primers,
            search_length,
            algorithm,
            edit_distance,
            max_mismatch,
            min_overlap,
            myers_patterns,
        )?;

        Ok(Self {
            matcher,
            local_output: String::new(),
            global_output: Arc::new(Mutex::new(output)),
        })
    }
}

impl<R: Record> ParallelProcessor<R> for PrimerTrimmer {
    fn process_record(&mut self, record: R) -> paraseq::parallel::Result<()> {
        // Convert sequence and quality to &str (keep bindings for lifetime)
        let seq_bytes = record.seq();
        let seq = std::str::from_utf8(&seq_bytes).map_err(|e| e.into_process_error())?;
        let qual_bytes = record
            .qual()
            .ok_or_else(|| anyhow::anyhow!("Missing quality scores"))?;
        let qual = std::str::from_utf8(qual_bytes).map_err(|e| e.into_process_error())?;

        // Search for paired primers and trim if found
        let search_result = self.matcher.search_primers(seq);

        if let Some((trimmed_seq, trimmed_qual)) = trim_sequence(seq, qual, &search_result) {
            // Write trimmed record to local buffer
            use std::fmt::Write;
            writeln!(self.local_output, "@{}", record.id_str())
                .map_err(|e| e.into_process_error())?;
            writeln!(self.local_output, "{}", trimmed_seq).map_err(|e| e.into_process_error())?;
            writeln!(self.local_output, "+").map_err(|e| e.into_process_error())?;
            writeln!(self.local_output, "{}", trimmed_qual).map_err(|e| e.into_process_error())?;
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        let mut global_out = self.global_output.lock().unwrap();
        global_out.write_all(self.local_output.as_bytes())?;
        global_out.flush()?;
        self.local_output.clear();
        Ok(())
    }
}
