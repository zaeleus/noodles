mod index;

use std::io::{self, Read};

use self::index::read_index;
use crate::bai::Index;

/// A BAM index (BAI) reader.
///
/// A BAM index has three top-level fields:
///
///   1. a magic number,
///   2. a list of reference sequences,
///   3. and optionally, the number of unmapped reads in the associated BAM.
///
/// While these fields can be read individually, consider using [`crate::bai::fs::read`] to read
/// the entire index at once.
///
/// # Examples
///
/// ```no_run
///# use std::{fs::File, io};
/// use noodles_bam::bai;
/// let mut reader = File::open("sample.bam.bai").map(bai::io::Reader::new)?;
/// let index = reader.read_index()?;
/// # Ok::<(), io::Error>(())
/// ```
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
    /// use noodles_bam::bai;
    /// let reader = bai::io::Reader::new(io::empty());
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
    /// use noodles_bam::bai;
    /// let mut reader = bai::io::Reader::new(io::empty());
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
    /// use noodles_bam::bai;
    /// let reader = bai::io::Reader::new(io::empty());
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
    /// Creates a BAM index reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_bam::bai;
    /// let reader = bai::io::Reader::new(io::empty());
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads the BAM index.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_bam::bai;
    /// let mut reader = File::open("sample.bam.bai").map(bai::io::Reader::new)?;
    /// let index = reader.read_index()?;
    /// # Ok::<(), std::io::Error>(())
    /// ```
    pub fn read_index(&mut self) -> io::Result<Index> {
        read_index(&mut self.inner)
    }
}
