use super::{Bin, Chunk};

#[derive(Debug, Default)]
pub struct Builder {
    id: u32,
    chunks: Vec<Chunk>,
}

impl Builder {
    pub fn set_id(&mut self, id: u32) -> &mut Self {
        self.id = id;
        self
    }

    pub fn add_chunk(&mut self, chunk: Chunk) -> &mut Self {
        if let Some(last_chunk) = self.chunks.last_mut() {
            if chunk.start() <= last_chunk.end() {
                *last_chunk = Chunk::new(last_chunk.start(), chunk.end());
                return self;
            }
        }

        self.chunks.push(chunk);

        self
    }

    pub fn build(self) -> Bin {
        Bin {
            id: self.id,
            chunks: self.chunks,
        }
    }
}
