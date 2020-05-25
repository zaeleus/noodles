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
}
