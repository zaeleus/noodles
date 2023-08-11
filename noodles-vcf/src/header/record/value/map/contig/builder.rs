use super::Contig;
use crate::header::record::value::map::{self, builder::BuildError};

#[derive(Default)]
pub struct Builder {
    length: Option<usize>,
    md5: Option<String>,
    url: Option<String>,
    idx: Option<usize>,
}

impl map::builder::Inner<Contig> for Builder {
    fn build(self) -> Result<Contig, BuildError> {
        Ok(Contig {
            length: self.length,
            md5: self.md5,
            url: self.url,
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

    /// Sets the MD5 hexdigest.
    pub fn set_md5(mut self, md5: String) -> Self {
        self.inner.md5 = Some(md5);
        self
    }

    /// Set the URL.
    pub fn set_url<U>(mut self, url: U) -> Self
    where
        U: Into<String>,
    {
        self.inner.url = Some(url.into());
        self
    }
}
