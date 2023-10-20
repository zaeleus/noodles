use std::{fs::File, io, path::Path};

use noodles_bgzf as bgzf;

use super::Reader;

/// A BCF reader builder.
pub struct Builder;

impl Builder {
    /// Builds a BCF reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_bcf::reader::Builder;
    /// let reader = Builder.build_from_path("sample.bcf")?;
    /// # Ok::<_, std::io::Error>(())
    /// ````
    pub fn build_from_path<P>(self, src: P) -> io::Result<Reader<bgzf::Reader<File>>>
    where
        P: AsRef<Path>,
    {
        File::open(src).map(Reader::new)
    }
}
