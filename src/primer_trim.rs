use crate::preparing_input::{ExpandedPrimerSet, MyersPatternSet, PrimerSet, reverse_complement};
use crate::primer_search::{PairedPrimerSearchResult, PrimerMatcher};
use crate::search_algos::Algorithm;
use paraseq::Record;
use paraseq::parallel::{IntoProcessError, ParallelProcessor};
use std::borrow::Cow;
use std::io::Write;
use std::sync::Arc;
use std::sync::Mutex;

/// Convert byte slice to string without UTF-8 validation.
///
/// SAFETY: This function is safe because:
/// - FASTQ sequences contain only A, C, G, T, N (all ASCII)
/// - FASTQ quality strings contain Phred+33 values (all ASCII)
/// - The paraseq parser validates UTF-8 during parsing before we see the data
/// - If parsing succeeded, the data is guaranteed to be valid UTF-8
#[inline]
fn bytes_to_str_unchecked(bytes: &[u8]) -> &str {
    unsafe { std::str::from_utf8_unchecked(bytes) }
}

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

    // Trim from both ends first (using original read positions)
    // trim_start: position to start keeping (trim everything before this)
    // trim_end: position to stop keeping (trim everything from this position to end)

    // Ensure we have valid trim positions
    let trim_end = result.trim_end.min(seq.len());
    if trim_end <= result.trim_start {
        // Invalid trimming positions or empty result
        return None;
    }

    // Trim the sequence and quality
    let trimmed_seq = &seq[result.trim_start..trim_end];
    let trimmed_qual = &qual[result.trim_start..trim_end];

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
pub struct AmpliconRecordProcessor {
    matcher: PrimerMatcher,
    local_output: String,
    global_output: Arc<Mutex<Box<dyn Write + Send>>>,
}

impl AmpliconRecordProcessor {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        output: Box<dyn Write + Send>,
        primers: PrimerSet,
        search_length: usize,
        algorithm: Algorithm,
        error_rate: f64,
        min_overlap: usize,
        myers_patterns: Option<MyersPatternSet>,
        expanded_primers: Option<ExpandedPrimerSet>,
    ) -> anyhow::Result<Self> {
        let matcher = PrimerMatcher::new(
            primers,
            search_length,
            algorithm,
            error_rate,
            min_overlap,
            myers_patterns,
            expanded_primers,
        )?;

        Ok(Self {
            matcher,
            local_output: String::new(),
            global_output: Arc::new(Mutex::new(output)),
        })
    }
}

impl<R: Record> ParallelProcessor<R> for AmpliconRecordProcessor {
    fn process_record(&mut self, record: R) -> paraseq::parallel::Result<()> {
        // Convert sequence and quality to &str using unchecked conversion
        // Safe: paraseq parser already validated UTF-8 during parsing
        let seq_bytes = record.seq();
        let seq = bytes_to_str_unchecked(&seq_bytes);
        let qual_bytes = record
            .qual()
            .ok_or_else(|| anyhow::anyhow!("Missing quality scores"))?;
        let qual = bytes_to_str_unchecked(qual_bytes);

        let search_result = self.matcher.search_primers(seq);

        let trimmed = trim_sequence(seq, qual, &search_result);

        if let Some((trimmed_seq, trimmed_qual)) = trimmed {
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
        let mut global_out = self
            .global_output
            .lock()
            .expect("Failed to acquire write lock on output");
        global_out.write_all(self.local_output.as_bytes())?;
        global_out.flush()?;
        self.local_output.clear();
        Ok(())
    }
}
