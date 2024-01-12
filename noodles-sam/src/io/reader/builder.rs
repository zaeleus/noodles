use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read},
    path::Path,
};

use noodles_bgzf as bgzf;

use super::Reader;
use crate::io::CompressionMethod;

/// A SAM reader builder.
#[derive(Debug, Default)]
pub struct Builder {
    compression_method: Option<CompressionMethod>,
}

impl Builder {
    /// Sets the compression method.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::io::{reader::Builder, CompressionMethod};
    /// let builder = Builder::default().set_compression_method(CompressionMethod::Bgzf);
    /// ```
    pub fn set_compression_method(mut self, compression_method: CompressionMethod) -> Self {
        self.compression_method = Some(compression_method);
        self
    }

    /// Builds a SAM reader from a path.
    ///
    /// By default, the compression method will be autodetected. This can be overridden by using
    /// [`Self::set_compression_method`].
    pub fn build_from_path<P>(mut self, src: P) -> io::Result<Reader<Box<dyn BufRead>>>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        if self.compression_method.is_none() {
            self.compression_method = match src.extension().and_then(|ext| ext.to_str()) {
                Some("gz" | "bgz") => Some(CompressionMethod::Bgzf),
                _ => Some(CompressionMethod::None),
            };
        }

        let reader = File::open(src)?;
        self.build_from_reader(reader)
    }

    /// Builds a SAM reader from a reader.
    pub fn build_from_reader<'r, R>(self, reader: R) -> io::Result<Reader<Box<dyn BufRead + 'r>>>
    where
        R: Read + 'r,
    {
        let inner: Box<dyn BufRead> = match self.compression_method {
            Some(CompressionMethod::Bgzf) => Box::new(bgzf::Reader::new(reader)),
            Some(CompressionMethod::None) | None => Box::new(BufReader::new(reader)),
        };

        Ok(Reader::new(inner))
    }
}
