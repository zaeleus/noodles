//! Alignment reader.

pub(crate) mod builder;
mod inner;

use std::io::{self, Read};

use noodles_sam as sam;

pub use self::builder::Builder;
use self::inner::Inner;

/// An alignment reader.
pub struct Reader<R>(Inner<R>);

impl<R> Reader<R>
where
    R: Read,
{
    /// Reads and parses an alignment header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment;
    /// let mut reader = alignment::io::reader::Builder::default().build_from_reader(io::empty())?;
    /// let header = reader.read_header()?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<sam::Header> {
        self.0.read_header()
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment;
    ///
    /// let mut reader = alignment::io::reader::Builder::default().build_from_reader(io::empty())?;
    /// let header = reader.read_header()?;
    /// let mut records = reader.records(&header);
    ///
    /// for result in records {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
    ) -> impl Iterator<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'r {
        self.0.records(header)
    }
}
