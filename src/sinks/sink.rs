pub trait RecordSink: Send + Clone + 'static {
    fn accept(&mut self, id: &str, seq: &str, qual: &str) -> std::io::Result<()>;
    fn on_batch_complete(&mut self) -> std::io::Result<()>;
}
