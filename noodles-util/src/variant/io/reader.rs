//! Variant reader.

pub(crate) mod builder;

pub use self::builder::Builder;

use std::io::{self, BufRead};

use noodles_vcf::{self as vcf, variant::RecordBuf};

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
    /// # use std::io::{self, Cursor};
    /// use noodles_vcf as vcf;
    /// use noodles_util::variant;
    ///
    /// let data = Cursor::new(b"##fileformat=VCFv4.4
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// sq0\t1\t.\tA\t.\t.\tPASS\t.
    /// ");
    ///
    /// let mut reader = variant::io::reader::Builder::default().build_from_reader(data)?;
    /// let actual = reader.read_header()?;
    ///
    /// let expected = vcf::Header::builder().build();
    ///
    /// assert_eq!(actual, expected);
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<vcf::Header> {
        self.inner.read_variant_header()
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Cursor};
    /// use noodles_vcf as vcf;
    /// use noodles_util::variant;
    ///
    /// let data = Cursor::new(b"##fileformat=VCFv4.3
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// sq0\t1\t.\tA\t.\t.\tPASS\t.
    /// ");
    ///
    /// let mut reader = variant::io::reader::Builder::default().build_from_reader(data)?;
    /// let header = reader.read_header()?;
    ///
    /// let mut records = reader.records(&header);
    ///
    /// assert!(records.next().transpose()?.is_some());
    /// assert!(records.next().is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records<'a>(
        &'a mut self,
        header: &'a vcf::Header,
    ) -> impl Iterator<Item = io::Result<RecordBuf>> + 'a {
        self.inner.variant_records(header)
    }
}
