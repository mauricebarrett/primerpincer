use super::RecordSink;

#[derive(Clone)]
pub struct SizeFilterSink<S: RecordSink> {
    inner: S,
    min_len: Option<usize>,
    max_len: Option<usize>,
}

impl<S: RecordSink> SizeFilterSink<S> {
    pub fn new(inner: S, min_len: Option<usize>, max_len: Option<usize>) -> Self {
        Self {
            inner,
            min_len,
            max_len,
        }
    }
}

impl<S: RecordSink> RecordSink for SizeFilterSink<S> {
    fn accept(&mut self, id: &str, seq: &str, qual: &str) -> std::io::Result<()> {
        let len = seq.len();
        if let Some(min) = self.min_len
            && len < min
        {
            return Ok(());
        }
        if let Some(max) = self.max_len
            && len > max
        {
            return Ok(());
        }
        self.inner.accept(id, seq, qual)
    }

    fn on_batch_complete(&mut self) -> std::io::Result<()> {
        self.inner.on_batch_complete()
    }
}
