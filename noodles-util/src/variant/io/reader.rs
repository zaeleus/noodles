//! Variant reader.

pub(crate) mod builder;

use std::io::{self, BufRead};

use noodles_bcf as bcf;
use noodles_vcf::{self as vcf, variant::io::Read};

pub use self::builder::Builder;
use crate::variant::Record;

pub(crate) enum Inner<R> {
    Vcf(vcf::io::Reader<R>),
    Bcf(bcf::io::Reader<R>),
}

/// A variant reader.
pub struct Reader<R> {
    pub(crate) inner: Inner<R>,
}

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
        match &mut self.inner {
            Inner::Vcf(reader) => reader.read_variant_header(),
            Inner::Bcf(reader) => reader.read_variant_header(),
        }
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
        match &mut self.inner {
            Inner::Vcf(reader) => {
                if !matches!(record, Record::Vcf(_)) {
                    *record = Record::Vcf(vcf::Record::default());
                }

                if let Record::Vcf(r) = record {
                    reader.read_record(r)
                } else {
                    unreachable!();
                }
            }
            Inner::Bcf(reader) => {
                if !matches!(record, Record::Bcf(_)) {
                    *record = Record::Bcf(bcf::Record::default());
                }

                if let Record::Bcf(r) = record {
                    reader.read_record(r)
                } else {
                    unreachable!();
                }
            }
        }
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
        match &mut self.inner {
            Inner::Vcf(reader) => reader.variant_records(header),
            Inner::Bcf(reader) => reader.variant_records(header),
        }
    }
}
