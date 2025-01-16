mod index;

use std::io::{self, Read};

use self::index::read_index;
use crate::gzi::Index;

/// A gzip index (GZI) reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf::gzi;
    /// let reader = gzi::io::Reader::new(io::empty());
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf::gzi;
    /// let mut reader = gzi::io::Reader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf::gzi;
    /// let reader = gzi::io::Reader::new(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a gzip index (GZI) reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::gzi;
    /// let data = [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];
    /// let reader = gzi::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a gzip index.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_bgzf::gzi;
    /// let mut reader = File::open("in.gzi").map(gzi::io::Reader::new)?;
    /// let index = reader.read_index()?;
    /// # Ok::<(), std::io::Error>(())
    /// ```
    pub fn read_index(&mut self) -> io::Result<Index> {
        read_index(&mut self.inner)
    }
}
