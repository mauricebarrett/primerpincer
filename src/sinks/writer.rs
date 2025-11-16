use super::RecordSink;
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

#[derive(Clone)]
pub struct WriterSink {
    local_output: String,
    global_output: Arc<Mutex<Box<dyn Write + Send>>>,
    written_count: Arc<AtomicUsize>,
}

impl WriterSink {
    pub fn new(output: Box<dyn Write + Send>, written_count: Arc<AtomicUsize>) -> Self {
        Self {
            local_output: String::new(),
            global_output: Arc::new(Mutex::new(output)),
            written_count,
        }
    }
}

impl RecordSink for WriterSink {
    fn accept(&mut self, id: &str, seq: &str, qual: &str) -> std::io::Result<()> {
        // Build FASTQ record into the local buffer without fallible formatting
        self.local_output.push('@');
        self.local_output.push_str(id);
        self.local_output.push('\n');
        self.local_output.push_str(seq);
        self.local_output.push('\n');
        self.local_output.push('+');
        self.local_output.push('\n');
        self.local_output.push_str(qual);
        self.local_output.push('\n');
        self.written_count.fetch_add(1, Ordering::Relaxed);
        Ok(())
    }

    fn on_batch_complete(&mut self) -> std::io::Result<()> {
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
