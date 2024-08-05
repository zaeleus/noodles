//! BED reader.

mod builder;
mod record;

pub use self::builder::Builder;
use std::io::{self, BufRead};

use self::record::{read_record_3, read_record_4, read_record_5, read_record_6};
use crate::Record;

/// A BED reader.
pub struct Reader<const N: usize, R> {
    inner: R,
}

impl<R, const N: usize> Reader<N, R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let data = [];
    /// let reader = bed::io::Reader::<3, _>::new(&data[..]);
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
    /// let mut reader = bed::io::Reader::<3, _>::new(&data[..]);
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
    /// let reader = bed::io::Reader::<3, _>::new(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<const N: usize, R> Reader<N, R>
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
    /// let reader = bed::io::Reader::<3, _>::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }
}

impl<R> Reader<3, R>
where
    R: BufRead,
{
    /// Reads a BED3+ record.
    pub fn read_record(&mut self, record: &mut Record<3>) -> io::Result<usize> {
        read_record_3(&mut self.inner, record)
    }
}

impl<R> Reader<4, R>
where
    R: BufRead,
{
    /// Reads a BED4+ record.
    pub fn read_record(&mut self, record: &mut Record<4>) -> io::Result<usize> {
        read_record_4(&mut self.inner, record)
    }
}

impl<R> Reader<5, R>
where
    R: BufRead,
{
    /// Reads a BED5+ record.
    pub fn read_record(&mut self, record: &mut Record<5>) -> io::Result<usize> {
        read_record_5(&mut self.inner, record)
    }
}

impl<R> Reader<6, R>
where
    R: BufRead,
{
    /// Reads a BED6+ record.
    pub fn read_record(&mut self, record: &mut Record<6>) -> io::Result<usize> {
        read_record_6(&mut self.inner, record)
    }
}
