pub(crate) mod index;

use std::io::{self, Read};

use noodles_bgzf as bgzf;

use self::index::read_index;
use super::Index;

/// A CSI reader.
pub struct Reader<R> {
    inner: bgzf::Reader<R>,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a CSI reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_csi as csi;
    /// let reader = File::open("sample.bcf.csi").map(csi::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner: bgzf::Reader::new(inner),
        }
    }

    /// Reads a CSI index.
    ///
    /// The position of the stream is expected to be at the beginning.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_csi as csi;
    /// let mut reader = File::open("sample.bcf.csi").map(csi::Reader::new)?;
    /// let index = reader.read_index();
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_index(&mut self) -> io::Result<Index> {
        read_index(&mut self.inner)
    }
}
