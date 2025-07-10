//! BCF writer.

mod builder;
pub(crate) mod header;
mod record;

use std::io::{self, Write};

use byteorder::WriteBytesExt;
use noodles_bgzf as bgzf;
use noodles_vcf::{self as vcf, header::StringMaps};

pub use self::builder::Builder;
use self::header::write_header;
pub(crate) use self::record::write_record;
use crate::Record;

pub(crate) const MAJOR: u8 = 2;
pub(crate) const MINOR: u8 = 2;

/// A BCF writer.
pub struct Writer<W> {
    inner: W,
    string_maps: StringMaps,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// let writer = bcf::io::Writer::from(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// let mut writer = bcf::io::Writer::from(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// let mut writer = bcf::io::Writer::from(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }
}

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
    /// use noodles_bcf as bcf;
    /// use noodles_vcf as vcf;
    ///
    /// let mut writer = bcf::io::Writer::new(io::sink());
    ///
    /// let header = vcf::Header::default();
    /// writer.write_header(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        write_file_format(&mut self.inner)?;

        self.string_maps = StringMaps::try_from(header)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        write_header(&mut self.inner, header)
    }

    /// Writes a record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// use noodles_core::Position;
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::{
    ///         record::value::{map::Contig, Map},
    ///         StringMaps,
    ///     },
    /// };
    ///
    /// let mut writer = bcf::io::Writer::new(io::sink());
    ///
    /// let mut header = vcf::Header::builder()
    ///     .add_contig("sq0", Map::<Contig>::new())
    ///     .build();
    /// *header.string_maps_mut() = StringMaps::try_from(&header)?;
    ///
    /// writer.write_header(&header)?;
    ///
    /// let record = bcf::Record::default();
    /// writer.write_record(&header, &record)?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn write_record(&mut self, header: &vcf::Header, record: &Record) -> io::Result<()> {
        write_record(&mut self.inner, header, &self.string_maps, record)
    }
}

impl<W> Writer<bgzf::io::Writer<W>>
where
    W: Write,
{
    /// Creates a BCF writer with a default compression level.
    ///
    /// The given stream is wrapped in a BGZF encoder.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// let writer = bcf::io::Writer::new(io::sink());
    /// ```
    pub fn new(writer: W) -> Self {
        Self::from(bgzf::io::Writer::new(writer))
    }

    /// Attempts to finish the output stream.
    ///
    /// This is typically only manually called if the underlying stream is needed before the writer
    /// is dropped.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// let mut writer = bcf::io::Writer::new(io::sink());
    /// writer.try_finish()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_finish(&mut self) -> io::Result<()> {
        self.inner.try_finish()
    }
}

impl<W> From<W> for Writer<W> {
    fn from(inner: W) -> Self {
        Self {
            inner,
            string_maps: StringMaps::default(),
        }
    }
}

impl<W> vcf::variant::io::Write for Writer<W>
where
    W: Write,
{
    fn write_variant_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        self.write_header(header)
    }

    fn write_variant_record(
        &mut self,
        header: &vcf::Header,
        record: &dyn vcf::variant::Record,
    ) -> io::Result<()> {
        write_record(&mut self.inner, header, &self.string_maps, record)
    }
}

fn write_file_format<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    use crate::io::MAGIC_NUMBER;

    writer.write_all(&MAGIC_NUMBER)?;
    writer.write_u8(MAJOR)?;
    writer.write_u8(MINOR)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_file_format() -> io::Result<()> {
        let mut buf = Vec::new();
        write_file_format(&mut buf)?;

        let expected = [
            b'B', b'C', b'F', // magic
            0x02, // major
            0x02, // minor
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
