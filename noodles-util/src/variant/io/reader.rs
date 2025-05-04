//! Variant reader.

pub(crate) mod builder;

pub use self::builder::Builder;

use std::io::{self, BufRead};

use noodles_bcf as bcf;
use noodles_vcf::{
    self as vcf,
    variant::{io::Read, Record},
};

enum Inner<R> {
    Vcf(vcf::io::Reader<R>),
    Bcf(bcf::io::Reader<R>),
}

/// A variant reader.
pub struct Reader<R> {
    inner: Inner<R>,
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
    ) -> impl Iterator<Item = io::Result<Box<dyn Record>>> + 'a {
        match &mut self.inner {
            Inner::Vcf(reader) => reader.variant_records(header),
            Inner::Bcf(reader) => reader.variant_records(header),
        }
    }
}
