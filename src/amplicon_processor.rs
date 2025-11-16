use crate::preparing_input::{ExpandedPrimerSet, MyersPatternSet, PrimerSet, reverse_complement};
use crate::primer_search::{PairedPrimerSearchResult, PrimerMatcher};
use crate::search_algos::Algorithm;
use crate::sinks::RecordSink;
use paraseq::Record;
use paraseq::parallel::{IntoProcessError, ParallelProcessor};
use std::borrow::Cow;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

#[inline]
fn bytes_to_str_unchecked(bytes: &[u8]) -> &str {
    unsafe { std::str::from_utf8_unchecked(bytes) }
}

fn trim_sequence<'a>(
    seq: &'a str,
    qual: &'a str,
    result: &PairedPrimerSearchResult,
) -> Option<(Cow<'a, str>, Cow<'a, str>)> {
    if !result.found {
        return None;
    }
    let trim_end = result.trim_end.min(seq.len());
    if trim_end <= result.trim_start {
        return None;
    }
    let trimmed_seq = &seq[result.trim_start..trim_end];
    let trimmed_qual = &qual[result.trim_start..trim_end];
    if result.needs_reverse_complement {
        let seq_rc = reverse_complement(trimmed_seq);
        let mut qual_bytes = trimmed_qual.as_bytes().to_vec();
        qual_bytes.reverse();
        let qual_rc = unsafe { String::from_utf8_unchecked(qual_bytes) };
        Some((Cow::Owned(seq_rc), Cow::Owned(qual_rc)))
    } else {
        Some((Cow::Borrowed(trimmed_seq), Cow::Borrowed(trimmed_qual)))
    }
}

#[derive(Clone)]
pub struct AmpliconRecordProcessor<S: RecordSink> {
    matcher: PrimerMatcher,
    sink: S,
    input_count: Arc<AtomicUsize>,
    trimmed_count: Arc<AtomicUsize>,
}

impl<S: RecordSink> AmpliconRecordProcessor<S> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        sink: S,
        primers: PrimerSet,
        search_length: usize,
        algorithm: Algorithm,
        error_rate: f64,
        min_overlap: usize,
        myers_patterns: Option<MyersPatternSet>,
        expanded_primers: Option<ExpandedPrimerSet>,
        input_count: Arc<AtomicUsize>,
        trimmed_count: Arc<AtomicUsize>,
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
            sink,
            input_count,
            trimmed_count,
        })
    }
}

impl<R: Record, S: RecordSink> ParallelProcessor<R> for AmpliconRecordProcessor<S> {
    fn process_record(&mut self, record: R) -> paraseq::parallel::Result<()> {
        self.input_count.fetch_add(1, Ordering::Relaxed);
        let seq_bytes = record.seq();
        let seq = bytes_to_str_unchecked(&seq_bytes);
        let qual_bytes = record.qual().ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "Missing quality scores")
        })?;
        let qual = bytes_to_str_unchecked(qual_bytes);

        let search_result = self.matcher.search_primers(seq);
        if let Some((trimmed_seq, trimmed_qual)) = trim_sequence(seq, qual, &search_result) {
            self.trimmed_count.fetch_add(1, Ordering::Relaxed);
            self.sink
                .accept(record.id_str(), &trimmed_seq, &trimmed_qual)
                .map_err(IntoProcessError::into_process_error)?;
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        self.sink
            .on_batch_complete()
            .map_err(IntoProcessError::into_process_error)
    }
}
