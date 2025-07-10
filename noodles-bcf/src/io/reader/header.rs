//! BCF header reader.

mod format_version;
pub(crate) mod magic_number;
pub mod vcf_header;

use std::io::{self, BufRead, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::{self as vcf, header::StringMaps};

use self::{format_version::read_format_version, magic_number::read_magic_number};
use crate::io::MAGIC_NUMBER;

/// A BCF header reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: Read,
{
    pub(super) fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads the magic number.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_bcf as bcf;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// let magic_number = header_reader.read_magic_number()?;
    /// assert_eq!(magic_number, *b"BCF");
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_magic_number(&mut self) -> io::Result<[u8; MAGIC_NUMBER.len()]> {
        read_magic_number(&mut self.inner)
    }

    /// Reads the format version.
    ///
    /// The format is returned as a major-minor version pair.
    ///
    /// The position of the stream is expected to be directly after the magic number.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_bcf as bcf;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number()?;
    ///
    /// let format_version = header_reader.read_format_version()?;
    /// assert_eq!(format_version, (2, 2));
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_format_version(&mut self) -> io::Result<(u8, u8)> {
        read_format_version(&mut self.inner)
    }

    /// Returns a VCF header reader.
    ///
    /// The caller is responsible of discarding any extra padding in the header text, e.g., using
    /// [`vcf_header::Reader::discard_to_end`].
    ///
    /// The position of the stream is expected to be directly after the format version.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{io::Read, fs::File};
    /// use noodles_bcf as bcf;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number()?;
    /// header_reader.read_format_version()?;
    ///
    /// let mut raw_vcf_header_reader = header_reader.raw_vcf_header_reader()?;
    ///
    /// let mut buf = Vec::new();
    /// raw_vcf_header_reader.read_to_end(&mut buf)?;
    ///
    /// raw_vcf_header_reader.discard_to_end()?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn raw_vcf_header_reader(&mut self) -> io::Result<vcf_header::Reader<&mut R>> {
        let len = self.inner.read_u32::<LittleEndian>().map(u64::from)?;
        Ok(vcf_header::Reader::new(&mut self.inner, len))
    }
}

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<vcf::Header>
where
    R: Read,
{
    let mut header_reader = Reader::new(reader);
    read_header_inner(&mut header_reader)
}

fn read_header_inner<R>(reader: &mut Reader<R>) -> io::Result<vcf::Header>
where
    R: Read,
{
    reader
        .read_magic_number()
        .and_then(magic_number::validate)?;

    reader.read_format_version()?;

    let mut raw_vcf_header_reader = reader.raw_vcf_header_reader()?;
    read_vcf_header(&mut raw_vcf_header_reader)
}

fn read_vcf_header<R>(reader: &mut vcf_header::Reader<R>) -> io::Result<vcf::Header>
where
    R: Read,
{
    let mut parser = vcf::header::Parser::default();
    let mut string_maps = StringMaps::default();

    let mut buf = Vec::new();

    while read_line(reader, &mut buf)? != 0 {
        let entry = parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        string_maps
            .insert_entry(&entry)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    reader.discard_to_end()?;

    let mut header = parser
        .finish()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *header.string_maps_mut() = string_maps;

    Ok(header)
}

fn read_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    dst.clear();

    match reader.read_until(LINE_FEED, dst)? {
        0 => Ok(0),
        n => {
            if dst.ends_with(&[LINE_FEED]) {
                dst.pop();

                if dst.ends_with(&[CARRIAGE_RETURN]) {
                    dst.pop();
                }
            }

            Ok(n)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_header() -> io::Result<()> {
        use vcf::header::FileFormat;

        const NUL: u8 = 0x00;

        let raw_header = b"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
";

        let mut data = vec![
            b'B', b'C', b'F', // magic
            0x02, 0x02, // major_version, minor_version
        ];
        data.extend(61u32.to_le_bytes()); // l_text
        data.extend(raw_header); // text
        data.push(NUL);

        let mut reader = &data[..];
        let actual = read_header(&mut reader)?;

        let expected = vcf::Header::builder()
            .set_file_format(FileFormat::new(4, 3))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
