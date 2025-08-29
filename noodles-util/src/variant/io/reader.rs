//! Variant reader.

pub(crate) mod builder;
mod inner;

use std::io::{self, BufRead};

use noodles_vcf as vcf;

pub use self::builder::Builder;
use self::inner::Inner;
use crate::variant::Record;

/// A variant reader.
pub struct Reader<R>(Inner<R>);

impl<R> Reader<R>
where
    R: BufRead,
{
    /// Reads and parses a variant header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::io::reader::Builder;
    ///
    /// let data = b"##fileformat=VCFv4.5
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// ";
    ///
    /// let mut reader = Builder::default().build_from_reader(&data[..])?;
    /// let _header = reader.read_header()?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<vcf::Header> {
        self.0.read_header()
    }

    /// Reads a variant record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{self, io::reader::Builder};
    ///
    /// let data = b"##fileformat=VCFv4.5
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// ";
    ///
    /// let mut reader = Builder::default().build_from_reader(&data[..])?;
    /// let header = reader.read_header()?;
    ///
    /// let mut record = variant::Record::default();
    ///
    /// while reader.read_record(&mut record)? != 0 {
    ///     // ...
    /// }
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        self.0.read_record(record)
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::io::reader::Builder;
    ///
    /// let data = b"##fileformat=VCFv4.5
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// ";
    ///
    /// let mut reader = Builder::default().build_from_reader(&data[..])?;
    /// let header = reader.read_header()?;
    ///
    /// let mut records = reader.records(&header);
    ///
    /// for result in records {
    ///     let _record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn records<'a>(
        &'a mut self,
        header: &'a vcf::Header,
    ) -> impl Iterator<Item = io::Result<Box<dyn vcf::variant::Record>>> + 'a {
        self.0.records(header)
    }
}

impl<'a> Reader<Box<dyn BufRead + 'a>> {
    /// Creates a variant reader.
    ///
    /// This attempts to autodetect the compression method and format of the input.
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::variant;
    /// let mut reader = variant::io::Reader::new(io::empty())?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn new<R>(reader: R) -> io::Result<Self>
    where
        R: BufRead + 'a,
    {
        Builder::default().build_from_reader(reader)
    }
}
