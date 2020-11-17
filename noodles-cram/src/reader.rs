pub mod block;
pub mod compression_header;
mod container;
mod encoding;
pub mod record;
mod records;
pub mod slice;

pub use self::records::Records;

use std::{
    io::{self, Read, Seek, SeekFrom},
    str,
};

use byteorder::{LittleEndian, ReadBytesExt};

use super::{file_definition::Version, Container, FileDefinition, MAGIC_NUMBER};

/// A CRAM reader.
///
/// The CRAM format is comprised of four main parts: 1) a file definition, 2) a file header, 3) a
/// list of data containers, and 4) a end-of-file (EOF) container.
pub struct Reader<R>
where
    R: Read,
{
    inner: R,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a CRAM reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    /// let mut reader = File::open("sample.bam").map(cram::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(reader: R) -> Self {
        Self { inner: reader }
    }

    /// Reads the CRAM file definition.
    ///
    /// The CRAM magic number is also checked.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    /// let mut reader = File::open("sample.cram").map(cram::Reader::new)?;
    /// let file_definition = reader.read_file_definition()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_file_definition(&mut self) -> io::Result<FileDefinition> {
        let magic = read_magic(&mut self.inner)?;

        if magic != MAGIC_NUMBER {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid CRAM header",
            ));
        }

        let format = read_format(&mut self.inner)?;
        let file_id = read_file_id(&mut self.inner)?;

        Ok(FileDefinition::new(format, file_id))
    }

    /// Reads the raw SAM header.
    ///
    /// The position of the stream is expected to be at the CRAM header container, i.e., directly
    /// after the file definition.
    ///
    /// This returns the raw SAM header as a [`String`]. It can subsequently be parsed as a
    /// [`sam::Header`].
    ///
    /// [`String`]: https://doc.rust-lang.org/std/string/struct.String.html
    /// [`sam::Header`]: ../../noodles_sam/header/struct.Header.html
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::Reader::new)?;
    /// reader.read_file_definition()?;
    ///
    /// let header = reader.read_file_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_file_header(&mut self) -> io::Result<String> {
        let header_container = self.read_container()?;

        if let Some(block) = header_container.blocks().first() {
            let data = block.decompressed_data()?;
            let mut reader = &data[..];

            let _header_len = reader.read_i32::<LittleEndian>()?;

            str::from_utf8(reader)
                .map(|s| s.into())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        } else {
            Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid header container: missing block for SAM header",
            ))
        }
    }

    fn read_container(&mut self) -> io::Result<Container> {
        let header = container::read_header(&mut self.inner)?;

        let blocks_len = header.block_count() as usize;
        let mut blocks = Vec::with_capacity(blocks_len);

        for _ in 0..blocks_len {
            let block = block::read_block(&mut self.inner)?;
            blocks.push(block);
        }

        Ok(Container::new(header, blocks))
    }

    /// Returns a iterator over records starting from the current stream position.
    ///
    /// The stream is expected to be at the start of a data container.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::Reader::new)?;
    /// reader.read_file_definition()?;
    /// reader.read_file_header()?;
    ///
    /// for result in reader.records() {
    ///     let record = result?;
    ///     println!("{:?}", record);
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self)
    }
}

impl<R> Reader<R>
where
    R: Read + Seek,
{
    /// Returns the current position of the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Cursor};
    /// use noodles_cram as cram;
    /// let data = Cursor::new(Vec::new());
    /// let mut reader = cram::Reader::new(data);
    /// let position = reader.position()?;
    /// assert_eq!(position, 0);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn position(&mut self) -> io::Result<u64> {
        self.inner.seek(SeekFrom::Current(0))
    }
}

fn read_magic<R>(reader: &mut R) -> io::Result<[u8; 4]>
where
    R: Read,
{
    let mut buf = [0; 4];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}

fn read_format<R>(reader: &mut R) -> io::Result<Version>
where
    R: Read,
{
    let mut buf = [0; 2];
    reader.read_exact(&mut buf)?;
    Ok(Version::new(buf[0], buf[1]))
}

fn read_file_id<R>(reader: &mut R) -> io::Result<[u8; 20]>
where
    R: Read,
{
    let mut buf = [0; 20];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_file_definition() -> Vec<u8> {
        let mut data = MAGIC_NUMBER.to_vec();

        let format = [0x03, 0x00];
        data.extend_from_slice(&format);

        let file_id = [
            0x00, 0x68, 0xac, 0xf3, 0x06, 0x4d, 0xaa, 0x1e, 0x29, 0xa4, 0xa0, 0x8c, 0x56, 0xee,
            0x91, 0x9b, 0x91, 0x04, 0x21, 0x1f,
        ];
        data.extend_from_slice(&file_id);

        data
    }

    #[test]
    fn test_read_file_definition() -> io::Result<()> {
        let data = build_file_definition();
        let mut reader = Reader::new(&data[..]);

        let actual = reader.read_file_definition()?;
        let expected = FileDefinition::new(
            Version::new(3, 0),
            [
                0x00, 0x68, 0xac, 0xf3, 0x06, 0x4d, 0xaa, 0x1e, 0x29, 0xa4, 0xa0, 0x8c, 0x56, 0xee,
                0x91, 0x9b, 0x91, 0x04, 0x21, 0x1f,
            ],
        );

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_file_definition_with_invalid_magic_number() {
        let data = b"BAM\x01";
        let mut reader = Reader::new(&data[..]);
        assert!(reader.read_file_definition().is_err());
    }
}
