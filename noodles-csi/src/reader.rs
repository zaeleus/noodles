//! CSI reader.

pub mod index;

use std::io::{self, Read};

use noodles_bgzf as bgzf;

use self::index::read_index;
use super::Index;

/// A CSI reader.
pub struct Reader<R> {
    inner: bgzf::Reader<R>,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    /// let reader = csi::Reader::new(io::empty());
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &bgzf::Reader<R> {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    /// let mut reader = csi::Reader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut bgzf::Reader<R> {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    /// let reader = csi::Reader::new(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> bgzf::Reader<R> {
        self.inner
    }
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
        read_index(&mut self.inner).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}
