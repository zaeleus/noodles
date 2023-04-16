mod record;
mod string_map;
mod value;

use std::{
    ffi::CString,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_bgzf as bgzf;
use noodles_vcf as vcf;

use self::record::write_record;
use super::header::StringMaps;

const MAJOR: u8 = 2;
const MINOR: u8 = 2;

/// A BCF writer.
pub struct Writer<W> {
    inner: W,
    string_maps: StringMaps,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let writer = bcf::Writer::from(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let mut writer = bcf::Writer::from(Vec::new());
    /// assert!(writer.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let mut writer = bcf::Writer::from(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }

    /// Writes a BCF file format.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// let mut writer = bcf::Writer::new(Vec::new());
    /// writer.write_file_format()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_file_format(&mut self) -> io::Result<()> {
        write_file_format(&mut self.inner)
    }

    /// Writes a VCF header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// use noodles_vcf as vcf;
    ///
    /// let mut writer = bcf::Writer::new(Vec::new());
    ///
    /// let header = vcf::Header::default();
    /// writer.write_header(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        self.string_maps = StringMaps::try_from(header)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        write_header(&mut self.inner, header)
    }

    /// Writes a record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::{self as bcf, header::StringMaps};
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::record::value::{map::Contig, Map},
    ///     record::Position,
    /// };
    ///
    /// let mut writer = bcf::Writer::new(Vec::new());
    ///
    /// let header = vcf::Header::builder()
    ///     .add_contig("sq0".parse()?, Map::<Contig>::new())
    ///     .build();
    ///
    /// writer.write_header(&header)?;
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::from(8))
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// writer.write_record(&header, &record)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write_record(&mut self, header: &vcf::Header, record: &vcf::Record) -> io::Result<()> {
        write_record(&mut self.inner, header, &self.string_maps, record)
    }
}

impl<W> Writer<bgzf::Writer<W>>
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
    /// use noodles_bcf as bcf;
    /// let writer = bcf::Writer::new(Vec::new());
    /// ```
    pub fn new(writer: W) -> Self {
        Self::from(bgzf::Writer::new(writer))
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
    /// let mut writer = bcf::Writer::new(Vec::new());
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

impl<W> vcf::VariantWriter for Writer<W>
where
    W: Write,
{
    fn write_variant_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        write_file_format(&mut self.inner)?;
        self.write_header(header)
    }

    fn write_variant_record(
        &mut self,
        header: &vcf::Header,
        record: &vcf::Record,
    ) -> io::Result<()> {
        write_record(&mut self.inner, header, &self.string_maps, record)
    }
}

fn write_file_format<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    use super::MAGIC_NUMBER;

    writer.write_all(MAGIC_NUMBER)?;
    writer.write_u8(MAJOR)?;
    writer.write_u8(MINOR)?;

    Ok(())
}

fn write_header<W>(writer: &mut W, header: &vcf::Header) -> io::Result<()>
where
    W: Write,
{
    let raw_header = header.to_string();
    let c_raw_header =
        CString::new(raw_header).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let text = c_raw_header.as_bytes_with_nul();
    let l_text =
        u32::try_from(text.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    writer.write_u32::<LittleEndian>(l_text)?;
    writer.write_all(text)?;

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

    #[test]
    fn test_write_header() -> io::Result<()> {
        let mut buf = Vec::new();
        let header = vcf::Header::default();
        write_header(&mut buf, &header)?;

        let mut expected = vec![
            0x3d, 0x00, 0x00, 0x00, // l_text = 61
        ];

        let text = b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\0";
        expected.extend_from_slice(text);

        assert_eq!(buf, expected);

        Ok(())
    }
}
