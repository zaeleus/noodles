mod record;

use std::io::{self, BufRead};

use self::record::{read_record_3, read_record_4, read_record_5, read_record_6};
use crate::Record;

/// A BED reader.
pub struct Reader<R, const N: usize> {
    inner: R,
}

impl<R, const N: usize> Reader<R, N> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let data = [];
    /// let reader = bed::io::Reader::<_, 3>::new(&data[..]);
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
    /// let mut reader = bed::io::Reader::<_, 3>::new(&data[..]);
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
    /// let reader = bed::io::Reader::<_, 3>::new(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R, const N: usize> Reader<R, N>
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
    /// let reader = bed::io::Reader::<_, 3>::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }
}

impl<R> Reader<R, 3>
where
    R: BufRead,
{
    /// Reads a BED3+ record.
    pub fn read_record(&mut self, record: &mut Record<3>) -> io::Result<usize> {
        read_record_3(&mut self.inner, record)
    }
}

impl<R> Reader<R, 4>
where
    R: BufRead,
{
    /// Reads a BED4+ record.
    pub fn read_record(&mut self, record: &mut Record<4>) -> io::Result<usize> {
        read_record_4(&mut self.inner, record)
    }
}

impl<R> Reader<R, 5>
where
    R: BufRead,
{
    /// Reads a BED5+ record.
    pub fn read_record(&mut self, record: &mut Record<5>) -> io::Result<usize> {
        read_record_5(&mut self.inner, record)
    }
}

impl<R> Reader<R, 6>
where
    R: BufRead,
{
    /// Reads a BED6+ record.
    pub fn read_record(&mut self, record: &mut Record<6>) -> io::Result<usize> {
        read_record_6(&mut self.inner, record)
    }
}
