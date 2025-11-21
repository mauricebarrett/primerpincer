use super::RecordSink;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Clone)]
pub struct SizeFilterSink<S: RecordSink> {
    inner: S,
    min_len: Option<usize>,
    max_len: Option<usize>,
    filtered_count: Arc<AtomicUsize>,
}

impl<S: RecordSink> SizeFilterSink<S> {
    pub fn new(
        inner: S,
        min_len: Option<usize>,
        max_len: Option<usize>,
        filtered_count: Arc<AtomicUsize>,
    ) -> Self {
        Self {
            inner,
            min_len,
            max_len,
            filtered_count,
        }
    }
}

impl<S: RecordSink> RecordSink for SizeFilterSink<S> {
    fn accept(&mut self, id: &str, seq: &str, qual: &str) -> std::io::Result<()> {
        let len = seq.len();
        if let Some(min) = self.min_len
            && len < min
        {
            self.filtered_count.fetch_add(1, Ordering::Relaxed);
            return Ok(());
        }
        if let Some(max) = self.max_len
            && len > max
        {
            self.filtered_count.fetch_add(1, Ordering::Relaxed);
            return Ok(());
        }
        self.inner.accept(id, seq, qual)
    }

    fn on_batch_complete(&mut self) -> std::io::Result<()> {
        self.inner.on_batch_complete()
    }
}
