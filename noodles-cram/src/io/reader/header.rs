//! CRAM header reader.

pub mod container;
mod file_id;
mod format_version;
pub(crate) mod magic_number;

use std::io::{self, BufRead, BufReader, Read};

use noodles_sam as sam;

use self::{
    file_id::read_file_id, format_version::read_format_version, magic_number::read_magic_number,
};
use crate::{FileDefinition, MAGIC_NUMBER, file_definition::Version};

/// A CRAM header reader.
pub struct Reader<R> {
    inner: R,
    version: Version,
}

impl<R> Reader<R>
where
    R: Read,
{
    pub(super) fn new(inner: R) -> Self {
        Self {
            inner,
            version: Version::default(),
        }
    }

    /// Reads the magic number.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_cram as cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// let magic_number = header_reader.read_magic_number()?;
    /// assert_eq!(magic_number, *b"CRAM");
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_magic_number(&mut self) -> io::Result<[u8; MAGIC_NUMBER.len()]> {
        read_magic_number(&mut self.inner)
    }

    /// Reads the format version.
    ///
    /// The position of the stream is expected to be directly after the magic number.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_cram as cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number()?;
    ///
    /// let _format_version = header_reader.read_format_version()?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_format_version(&mut self) -> io::Result<Version> {
        let version = read_format_version(&mut self.inner)?;
        self.version = version;
        Ok(version)
    }

    /// Reads the file ID.
    ///
    /// The position of the stream is expected to be directly after the format version.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_cram as cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number()?;
    /// header_reader.read_format_version()?;
    ///
    /// let _file_id = header_reader.read_file_id()?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_file_id(&mut self) -> io::Result<[u8; 20]> {
        read_file_id(&mut self.inner)
    }

    /// Returns a header container reader.
    ///
    /// The caller is responsible of discarding any extra padding in the header container, e.g.,
    /// using [`container::Reader::discard_to_end`].
    ///
    /// The position of the stream is expected to be at the start of a header container.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_cram as cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number()?;
    /// header_reader.read_format_version()?;
    /// header_reader.read_file_id()?;
    ///
    /// let mut container_reader = header_reader.container_reader()?;
    /// let header = { /* ... */ };
    /// container_reader.discard_to_end()?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn container_reader(&mut self) -> io::Result<container::Reader<&mut R>> {
        let len = container::read_header(&mut self.inner, self.version)?;
        Ok(container::Reader::new(&mut self.inner, len, self.version))
    }
}

pub(super) fn read_file_definition<R>(reader: &mut R) -> io::Result<FileDefinition>
where
    R: Read,
{
    let mut header_reader = Reader::new(reader);
    read_file_definition_inner(&mut header_reader)
}

fn read_file_definition_inner<R>(reader: &mut Reader<R>) -> io::Result<FileDefinition>
where
    R: Read,
{
    reader
        .read_magic_number()
        .and_then(magic_number::validate)?;

    let version = reader.read_format_version()?;
    let file_id = reader.read_file_id()?;

    Ok(FileDefinition::new(version, file_id))
}

pub(super) fn read_file_header<R>(reader: &mut R, version: Version) -> io::Result<sam::Header>
where
    R: Read,
{
    let mut header_reader = Reader::new(reader);
    header_reader.version = version;
    read_file_header_inner(&mut header_reader)
}

fn read_file_header_inner<R>(reader: &mut Reader<R>) -> io::Result<sam::Header>
where
    R: Read,
{
    let mut container_reader = reader.container_reader()?;

    let header = {
        let mut raw_sam_header_reader = container_reader.raw_sam_header_reader()?;
        let header = read_sam_header(&mut raw_sam_header_reader)?;
        raw_sam_header_reader.discard_to_end()?;
        header
    };

    container_reader.discard_to_end()?;

    Ok(header)
}

fn read_sam_header<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: Read,
{
    let mut parser = sam::header::Parser::default();

    let mut header_reader = BufReader::new(reader);
    let mut buf = Vec::new();

    while read_line(&mut header_reader, &mut buf)? != 0 {
        parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    Ok(parser.finish())
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
    use crate::file_definition::Version;

    #[test]
    fn test_read_file_definition() -> Result<(), Box<dyn std::error::Error>> {
        let src = [
            0x43, 0x52, 0x41, 0x4d, // magic number = b"CRAM"
            0x03, 0x00, // format version = (3, 0)
            0x00, 0x68, 0xac, 0xf3, 0x06, 0x4d, 0xaa, 0x1e, 0x29, 0xa4, 0xa0, 0x8c, 0x56, 0xee,
            0x91, 0x9b, 0x91, 0x04, 0x21, 0x1f, // file ID
        ];

        let mut reader = &src[..];
        let actual = read_file_definition(&mut reader)?;

        let (_, file_id) = src.split_last_chunk().unwrap();
        let expected = FileDefinition::new(Version::new(3, 0), *file_id);

        assert_eq!(actual, expected);

        Ok(())
    }
}
