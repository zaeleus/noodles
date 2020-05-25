#[derive(Debug)]
pub struct Chunk {
    start: u64,
    end: u64,
}

impl Chunk {
    pub fn new(start: u64, end: u64) -> Self {
        Self { start, end }
    }
}
