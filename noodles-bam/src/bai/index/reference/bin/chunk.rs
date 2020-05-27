use noodles_bgzf as bgzf;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Chunk {
    chunk_beg: bgzf::VirtualPosition,
    chunk_end: bgzf::VirtualPosition,
}

impl Chunk {
    pub fn new(start: bgzf::VirtualPosition, end: bgzf::VirtualPosition) -> Self {
        Self {
            chunk_beg: start,
            chunk_end: end,
        }
    }

    pub fn start(&self) -> bgzf::VirtualPosition {
        self.chunk_beg
    }

    pub fn start_mut(&mut self) -> &mut bgzf::VirtualPosition {
        &mut self.chunk_beg
    }

    pub fn end(&self) -> bgzf::VirtualPosition {
        self.chunk_end
    }

    pub fn end_mut(&mut self) -> &mut bgzf::VirtualPosition {
        &mut self.chunk_end
    }
}
