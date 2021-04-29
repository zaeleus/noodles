pub mod record;
mod value;

use std::{
    convert::TryFrom,
    ffi::CStr,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf as bgzf;

use super::Record;

static MAGIC_NUMBER: &[u8] = b"BCF";

/// A BCF reader.
///
/// The BCF format is comprised of two parts: 1) a VCF header and 2) a list of records.
pub struct Reader<R> {
    inner: bgzf::Reader<R>,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a BCF reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    /// let mut reader = File::open("sample.bcf").map(bcf::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(reader: R) -> Self {
        Self {
            inner: bgzf::Reader::new(reader),
        }
    }

    /// Reads the BCF file format.
    ///
    /// The BCF magic number is also checked.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// This returns the major and minor format versions as a tuple.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    /// let mut reader = File::open("sample.bcf").map(bcf::Reader::new)?;
    /// let (major, minor) = reader.read_file_format()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_file_format(&mut self) -> io::Result<(u8, u8)> {
        let mut buf = [0; 5];

        self.inner.read_exact(&mut buf)?;

        let magic = &buf[..3];

        if magic != MAGIC_NUMBER {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid file format",
            ));
        }

        let major = buf[3];
        let minor = buf[4];

        Ok((major, minor))
    }

    /// Reads the raw VCF header.
    ///
    /// The position of the stream is expected to be directly after the file format.
    ///
    /// This returns the raw VCF header as a [`std::string::String`]. It can subsequently be parsed
    /// as a [`noodles_vcf::Header`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    /// let mut reader = File::open("sample.bcf").map(bcf::Reader::new)?;
    /// reader.read_file_format()?;
    /// let header = reader.read_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<String> {
        // § 6.2 Header (2021-01-13) doesn't actually describe how the header is stored. There
        // seems to be a 32-bit integer for the header length, which adds +1 for a null terminator.
        // This means the header string itself is null-terminated.

        let header_len = self.inner.read_i32::<LittleEndian>().and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let mut buf = vec![0; header_len];
        self.inner.read_exact(&mut buf)?;

        CStr::from_bytes_with_nul(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|c_header| {
                c_header
                    .to_str()
                    .map(|s| s.into())
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    /// Reads a single record.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// If successful, the record size is returned. If a record size of 0 is returned, the stream
    /// reached EOF.
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        let l_shared = match self.inner.read_u32::<LittleEndian>() {
            Ok(len) => len,
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
            Err(e) => return Err(e),
        };

        let l_indiv = self.inner.read_u32::<LittleEndian>()?;

        let record_len = l_shared
            .checked_add(l_indiv)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "invalid record length: l_shared = {}, l_indiv = {}",
                        l_shared, l_indiv
                    ),
                )
            })
            .and_then(|len| {
                usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })?;

        record.resize(record_len);
        self.inner.read_exact(record)?;

        Ok(record_len)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    fn compress(data: &[u8]) -> io::Result<Vec<u8>> {
        let mut writer = bgzf::Writer::new(Vec::new());
        writer.write_all(data)?;
        writer.finish()
    }

    #[test]
    fn test_read_file_format() -> io::Result<()> {
        let data = compress(b"BCF\x02\x01")?;
        let mut reader = Reader::new(&data[..]);

        let (major, minor) = reader.read_file_format()?;

        assert_eq!(major, 2);
        assert_eq!(minor, 1);

        Ok(())
    }

    #[test]
    fn test_read_file_format_with_an_invalid_magic_number() -> io::Result<()> {
        let data = compress(b"BAM\x02\x01")?;
        let mut reader = Reader::new(&data[..]);
        assert!(reader.read_file_format().is_err());
        Ok(())
    }
}
