mod record;
mod string_map;
mod value;

use std::{
    convert::TryFrom,
    ffi::CString,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_bgzf as bgzf;
use noodles_vcf as vcf;

use super::{header::StringMap, MAGIC_NUMBER};

const MAJOR: u8 = 2;
const MINOR: u8 = 2;

/// A BCF writer.
pub struct Writer<W>
where
    W: Write,
{
    inner: bgzf::Writer<W>,
}

impl<W> Writer<W>
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
        Self {
            inner: bgzf::Writer::new(writer),
        }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let writer = bcf::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        self.inner.get_ref()
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
        self.inner.write_all(MAGIC_NUMBER)?;
        self.inner.write_u8(MAJOR)?;
        self.inner.write_u8(MINOR)?;
        Ok(())
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
        let raw_header = header.to_string();
        let c_raw_header =
            CString::new(raw_header).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let text = c_raw_header.as_bytes_with_nul();
        let l_text = u32::try_from(text.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        self.inner.write_u32::<LittleEndian>(l_text)?;
        self.inner.write_all(text)?;

        Ok(())
    }

    /// Writes a VCF record.
    pub fn write_vcf_record(
        &mut self,
        header: &vcf::Header,
        string_map: &StringMap,
        record: &vcf::Record,
    ) -> io::Result<()> {
        let mut site_buf = Vec::new();
        record::write_site(&mut site_buf, header, string_map, record)?;

        let l_shared = u32::try_from(site_buf.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let mut genotypes_buf = Vec::new();

        if let Some(format) = record.format() {
            record::write_genotypes(&mut genotypes_buf, string_map, format, record.genotypes())?;
        };

        let l_indiv = u32::try_from(genotypes_buf.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        self.inner.write_u32::<LittleEndian>(l_shared)?;
        self.inner.write_u32::<LittleEndian>(l_indiv)?;
        self.inner.write_all(&site_buf)?;
        self.inner.write_all(&genotypes_buf)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::Reader;

    use super::*;

    #[test]
    fn test_write_file_format() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());
        writer.write_file_format()?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().as_slice());
        let file_format = reader.read_file_format()?;

        assert_eq!(file_format, (2, 2));

        Ok(())
    }

    #[test]
    fn test_write_header() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());
        writer.write_file_format()?;

        let header = vcf::Header::default();
        writer.write_header(&header)?;

        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().as_slice());
        reader.read_file_format()?;
        let actual = reader.read_header()?;

        let expected = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

        assert_eq!(actual, expected);

        Ok(())
    }
}
