use super::Contig;
use crate::header::record::value::map::{self, builder::BuildError};

#[derive(Default)]
pub struct Builder {
    length: Option<usize>,
    idx: Option<usize>,
}

impl map::builder::Inner<Contig> for Builder {
    fn build(self) -> Result<Contig, BuildError> {
        Ok(Contig {
            length: self.length,
            idx: self.idx,
        })
    }
}

impl map::builder::Indexed<Contig> for Builder {
    fn set_idx(mut self, idx: usize) -> Self {
        self.idx = Some(idx);
        self
    }
}

impl map::Builder<Contig> {
    /// Sets the length.
    pub fn set_length(mut self, length: usize) -> Self {
        self.inner.length = Some(length);
        self
    }
}
