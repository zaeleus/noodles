//! Variant writer.

mod builder;

pub use self::builder::Builder;

use std::io;

use noodles_vcf as vcf;

/// A variant writer.
pub struct Writer {
    inner: Box<dyn vcf::VariantWriter>,
}

impl Writer {
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
        self.inner.write_variant_header(header)
    }

    /// Writes a variant record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_vcf::{self as vcf, Record};
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
    /// let s = "sq0\t8\t.\tA\t.\t.\tPASS\t.";
    /// let record = Record::try_from_str(s, &header)?;
    /// writer.write_record(&header, &record)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write_record(&mut self, header: &vcf::Header, record: &vcf::Record) -> io::Result<()> {
        self.inner.write_variant_record(header, record)
    }
}
