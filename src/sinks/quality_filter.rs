use super::RecordSink;
use once_cell::sync::Lazy;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Convert Phred+33 quality score ASCII to numeric Q score
#[inline]
fn phred33_to_q(byte: u8) -> u8 {
    byte.saturating_sub(33)
}

const PHRED_TABLE_SIZE: usize = 94; // Support Q scores 0..93 (covers Phred+33 printable ASCII)

/// Lookup table for error probabilities of each possible Phred score.
/// This avoids calling powf for every base in a read.
static ERROR_PROB_TABLE: Lazy<[f64; PHRED_TABLE_SIZE]> = Lazy::new(|| {
    let mut table = [0.0_f64; PHRED_TABLE_SIZE];
    let mut q = 0;
    while q < PHRED_TABLE_SIZE {
        table[q] = 10f64.powf(-(q as f64) / 10.0);
        q += 1;
    }
    table
});

/// Calculate mean Phred score by:
/// 1. Convert each Phred score to error probability
/// 2. Average the error probabilities
/// 3. Convert back to Phred score
fn mean_phred_score(qual: &str) -> f64 {
    if qual.is_empty() {
        return 0.0;
    }

    // Step 1 & 2: Convert each Phred to error prob via lookup table, then sum
    let sum_error_prob: f64 = qual
        .bytes()
        .map(|b| {
            let q = phred33_to_q(b) as usize;
            let idx = q.min(PHRED_TABLE_SIZE - 1);
            ERROR_PROB_TABLE[idx]
        })
        .sum();

    // Average error probability
    let mean_error_prob = sum_error_prob / (qual.len() as f64);

    // Step 3: Convert back to Phred score
    if mean_error_prob <= 0.0 || mean_error_prob >= 1.0 {
        0.0 // Handle edge cases
    } else {
        -10.0 * mean_error_prob.log10()
    }
}

#[derive(Clone)]
pub struct QualityFilterSink<S: RecordSink> {
    inner: S,
    min_average_quality: Option<f64>,
    filtered_count: Arc<AtomicUsize>,
}

impl<S: RecordSink> QualityFilterSink<S> {
    pub fn new(
        inner: S,
        min_average_quality: Option<f64>,
        filtered_count: Arc<AtomicUsize>,
    ) -> Self {
        Self {
            inner,
            min_average_quality,
            filtered_count,
        }
    }
}

impl<S: RecordSink> RecordSink for QualityFilterSink<S> {
    fn accept(&mut self, id: &str, seq: &str, qual: &str) -> std::io::Result<()> {
        // If quality filtering is enabled, check the mean Phred score
        if let Some(threshold) = self.min_average_quality {
            let mean_q = mean_phred_score(qual);
            if mean_q < threshold {
                self.filtered_count.fetch_add(1, Ordering::Relaxed);
                return Ok(()); // Filter out - quality too low
            }
        }
        // Pass through to next sink
        self.inner.accept(id, seq, qual)
    }

    fn on_batch_complete(&mut self) -> std::io::Result<()> {
        self.inner.on_batch_complete()
    }
}
