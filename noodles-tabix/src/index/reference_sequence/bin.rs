mod chunk;

pub use self::chunk::Chunk;

#[derive(Debug)]
pub struct Bin {
    bin: u32,
    chunks: Vec<Chunk>,
}

impl Bin {
    pub fn new(bin: u32, chunks: Vec<Chunk>) -> Self {
        Self { bin, chunks }
    }

    pub fn bin(&self) -> u32 {
        self.bin
    }

    pub fn chunks(&self) -> &[Chunk] {
        &self.chunks
    }
}
