//! Variant reader.

pub(crate) mod builder;

pub use self::builder::Builder;

use std::io::{self, BufRead};

use noodles_vcf::{self as vcf, variant::Record};

/// A variant reader.
pub struct Reader<R> {
    inner: Box<dyn vcf::variant::io::Read<R>>,
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
        self.inner.read_variant_header()
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
        self.inner.variant_records(header)
    }
}
