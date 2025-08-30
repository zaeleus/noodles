//! Variant writer.

pub(crate) mod builder;
mod inner;

use std::io::{self, Write};

use noodles_vcf as vcf;

pub use self::builder::Builder;
use self::inner::Inner;

/// A variant writer.
pub struct Writer<W>(Inner<W>)
where
    W: Write;

impl<W> Writer<W>
where
    W: Write,
{
    /// Writes a VCF header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_vcf as vcf;
    /// use noodles_util::variant::{self, io::{CompressionMethod, Format}};
    ///
    /// let mut writer = variant::io::writer::Builder::default()
    ///     .set_format(Format::Bcf)
    ///     .set_compression_method(Some(CompressionMethod::Bgzf))
    ///     .build_from_writer(io::sink());
    ///
    /// let header = vcf::Header::default();
    /// writer.write_header(&header)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        self.0.write_header(header)
    }

    /// Writes a variant record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_vcf as vcf;
    /// use noodles_util::variant::{self, io::Format};
    ///
    /// let mut writer = variant::io::writer::Builder::default()
    ///     .set_format(Format::Vcf)
    ///     .set_compression_method(None)
    ///     .build_from_writer(io::sink());
    ///
    /// let header = vcf::Header::default();
    /// writer.write_header(&header)?;
    ///
    /// let record = vcf::variant::RecordBuf::default();
    /// writer.write_record(&header, &record)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write_record(
        &mut self,
        header: &vcf::Header,
        record: &dyn vcf::variant::Record,
    ) -> io::Result<()> {
        self.0.write_record(header, record)
    }
}
