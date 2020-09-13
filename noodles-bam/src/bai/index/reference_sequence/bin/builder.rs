use super::{Bin, Chunk};

#[derive(Debug, Default)]
pub struct Builder {
    bin: u32,
    chunks: Vec<Chunk>,
}

impl Builder {
    pub fn set_bin(&mut self, bin: u32) -> &mut Self {
        self.bin = bin;
        self
    }

    pub fn add_chunk(&mut self, chunk: Chunk) -> &mut Self {
        if let Some(last_chunk) = self.chunks.last_mut() {
            if chunk.start() <= last_chunk.end() {
                *last_chunk = Chunk::new(last_chunk.start(), chunk.end());
            }
        } else {
            self.chunks.push(chunk);
        }

        self
    }

    pub fn build(self) -> Bin {
        Bin {
            bin: self.bin,
            chunks: self.chunks,
        }
    }
}
