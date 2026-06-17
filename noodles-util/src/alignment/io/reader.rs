//! Alignment reader.

pub(crate) mod builder;
mod inner;

use std::io::{self, Read};

use noodles_sam as sam;

pub use self::builder::Builder;
use self::inner::Inner;
use crate::alignment::Record;

/// An alignment reader.
pub struct Reader<R>(Inner<R>);

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates an alignment reader.
    ///
    /// This attempts to autodetect the compression method and format of the input.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment;
    /// let reader = alignment::io::Reader::new(io::empty())?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn new(reader: R) -> io::Result<Self> {
        Builder::default().build_from_reader(reader)
    }

    /// Reads and parses an alignment header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment;
    /// let mut reader = alignment::io::Reader::new(io::empty())?;
    /// let header = reader.read_header()?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<sam::Header> {
        self.0.read_header()
    }

    /// Reads an alignment record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment;
    ///
    /// let mut reader = alignment::io::Reader::new(io::empty())?;
    /// let header = reader.read_header()?;
    ///
    /// let mut record = alignment::Record::default();
    ///
    /// while reader.read_record(&header, &mut record)? != 0 {
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_record(&mut self, header: &sam::Header, record: &mut Record) -> io::Result<usize> {
        self.0.read_record(header, record)
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment;
    ///
    /// let mut reader = alignment::io::Reader::new(io::empty())?;
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
