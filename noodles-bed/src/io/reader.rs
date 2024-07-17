mod record;

use std::io::{self, BufRead};

use self::record::read_record;
use crate::Record;

/// A BED reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let data = [];
    /// let reader = bed::io::Reader::new(&data[..]);
    /// assert!(reader.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let data = [];
    /// let mut reader = bed::io::Reader::new(&data[..]);
    /// assert!(reader.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let data = [];
    /// let reader = bed::io::Reader::new(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Reader<R>
where
    R: BufRead,
{
    /// Creates a BED reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let data = [];
    /// let reader = bed::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a record.
    pub fn read_record<const N: usize>(&mut self, record: &mut Record<N>) -> io::Result<usize> {
        read_record(&mut self.inner, record)
    }
}
