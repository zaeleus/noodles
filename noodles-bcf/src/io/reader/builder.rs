use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

use noodles_bgzf as bgzf;

use super::Reader;
use crate::io::CompressionMethod;

/// A BCF reader builder.
#[derive(Default)]
pub struct Builder {
    compression_method: Option<CompressionMethod>,
}

impl Builder {
    /// Sets the compression method.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::io::{reader::Builder, CompressionMethod};
    /// let builder = Builder::default().set_compression_method(CompressionMethod::Bgzf);
    /// ```
    pub fn set_compression_method(mut self, compression_method: CompressionMethod) -> Self {
        self.compression_method = Some(compression_method);
        self
    }

    /// Builds a BCF reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_bcf::io::reader::Builder;
    /// let reader = Builder::default().build_from_path("sample.bcf")?;
    /// # Ok::<_, std::io::Error>(())
    /// ````
    pub fn build_from_path<P>(self, src: P) -> io::Result<Reader<Box<dyn Read>>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(src)?;
        self.build_from_reader(file)
    }

    /// Builds a BCF reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf::io::reader::Builder;
    /// let reader = Builder::default().build_from_reader(io::empty());
    /// ```
    pub fn build_from_reader<'r, R>(self, reader: R) -> io::Result<Reader<Box<dyn Read + 'r>>>
    where
        R: Read + 'r,
    {
        let inner: Box<dyn Read> = match self.compression_method {
            Some(CompressionMethod::Bgzf) | None => Box::new(bgzf::io::Reader::new(reader)),
            Some(CompressionMethod::None) => Box::new(reader),
        };

        Ok(Reader::from(inner))
    }
}
