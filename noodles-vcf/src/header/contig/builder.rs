use indexmap::IndexMap;

use super::{tag, Contig, Name, ID};

#[derive(Default)]
pub struct Builder {
    id: Option<Name>,
    len: Option<usize>,
    idx: Option<usize>,
    other_fields: IndexMap<tag::Other, String>,
}

impl Builder {
    pub fn set_id(mut self, id: Name) -> Self {
        self.id = Some(id);
        self
    }

    pub fn set_len(mut self, len: usize) -> Self {
        self.len = Some(len);
        self
    }

    pub fn set_idx(mut self, idx: usize) -> Self {
        self.idx = Some(idx);
        self
    }

    pub fn insert(mut self, key: tag::Other, value: String) -> Self {
        self.other_fields.insert(key, value);
        self
    }

    pub fn build(self) -> Result<Contig, BuildError> {
        Ok(Contig {
            id: self.id.ok_or(BuildError::MissingField(ID))?,
            len: self.len,
            idx: self.idx,
            fields: self.other_fields,
        })
    }
}

pub enum BuildError {
    MissingField(&'static str),
}
