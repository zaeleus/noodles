//! CRAM reader and record iterator.

mod builder;
pub(crate) mod container;
pub(crate) mod data_container;
pub(crate) mod header_container;
pub(crate) mod num;
mod query;
pub(crate) mod record;
mod records;

pub use self::{builder::Builder, query::Query, records::Records};

use std::io::{self, Read, Seek, SeekFrom};

use bytes::BytesMut;
use noodles_core::Region;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::{crai, file_definition::Version, FileDefinition, MAGIC_NUMBER};
use crate::data_container::DataContainer;

/// A CRAM reader.
///
/// The CRAM format is comprised of four main parts: 1) a file definition, 2) a file header, 3) a
/// list of data containers, and 4) an end-of-file (EOF) container.
///
/// # Examples
///
/// ```no_run
/// # use std::{fs::File, io};
/// use noodles_cram as cram;
/// use noodles_fasta as fasta;
///
/// let mut reader = File::open("sample.cram").map(cram::Reader::new)?;
/// reader.read_file_definition()?;
///
/// let header = reader.read_file_header()?;
///
/// for result in reader.records(&header) {
///     let record = result?;
///     println!("{:?}", record);
/// }
///
/// # Ok::<_, io::Error>(())
/// ```
pub struct Reader<R> {
    inner: R,
    reference_sequence_repository: fasta::Repository,
    buf: BytesMut,
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
    pub fn new(inner: R) -> Self {
        Builder::default().build_from_reader(inner)
    }

    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let data = [];
    /// let reader = cram::Reader::new(&data[..]);
    /// assert!(reader.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let data = [];
    /// let mut reader = cram::Reader::new(&data[..]);
    /// assert!(reader.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Unwraps and returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let data = [];
    /// let reader = cram::Reader::new(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }

    pub(crate) fn reference_sequence_repository(&self) -> &fasta::Repository {
        &self.reference_sequence_repository
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
        read_magic_number(&mut self.inner)?;

        let format = read_format(&mut self.inner)?;
        let file_id = read_file_id(&mut self.inner)?;

        Ok(FileDefinition::new(format, file_id))
    }

    /// Reads the raw SAM header.
    ///
    /// The position of the stream is expected to be at the CRAM header container, i.e., directly
    /// after the file definition.
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
    pub fn read_file_header(&mut self) -> io::Result<sam::Header> {
        use self::header_container::read_header_container;
        read_header_container(&mut self.inner, &mut self.buf)
    }

    pub(crate) fn read_data_container_with_container_header(
        &mut self,
    ) -> io::Result<Option<(crate::data_container::Header, DataContainer)>> {
        use self::data_container::read_data_container_with_container_header;
        read_data_container_with_container_header(&mut self.inner, &mut self.buf)
    }

    /// Reads a data container.
    ///
    /// This returns `None` if the container header is the EOF container header, which signals the
    /// end of the stream.
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
    /// while let Some(container) = reader.read_data_container()? {
    ///     // ...
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_data_container(&mut self) -> io::Result<Option<DataContainer>> {
        use self::data_container::read_data_container;

        read_data_container(&mut self.inner, &mut self.buf)
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
    /// use noodles_fasta as fasta;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::Reader::new)?;
    /// reader.read_file_definition()?;
    ///
    /// let header = reader.read_file_header()?;
    ///
    /// for result in reader.records(&header) {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records<'r>(&'r mut self, header: &'r sam::Header) -> Records<'r, R> {
        Records::new(self, header)
    }
}

impl<R> Reader<R>
where
    R: Read + Seek,
{
    /// Seeks the underlying reader to the given position.
    ///
    /// Positions typically come from the associated CRAM index file.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use std::io::SeekFrom;
    /// use noodles_cram as cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::Reader::new)?;
    /// reader.seek(SeekFrom::Start(17711))?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.inner.seek(pos)
    }

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
        self.inner.stream_position()
    }

    /// Returns an iterator over records that intersects the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram::{self as cram, crai};
    /// use noodles_fasta as fasta;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::Reader::new)?;
    /// reader.read_file_definition()?;
    ///
    /// let header = reader.read_file_header()?;
    /// let index = crai::read("sample.cram.crai")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<'a>(
        &'a mut self,
        header: &'a sam::Header,
        index: &'a crai::Index,
        region: &Region,
    ) -> io::Result<Query<'_, R>> {
        let reference_sequence_id = header
            .reference_sequences()
            .get_index_of(region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "invalid reference sequence name",
                )
            })?;

        Ok(Query::new(
            self,
            header,
            index,
            reference_sequence_id,
            region.interval(),
        ))
    }
}

impl<R> sam::AlignmentReader<R> for Reader<R>
where
    R: Read,
{
    fn read_alignment_header(&mut self) -> io::Result<sam::Header> {
        self.read_file_definition()?;
        self.read_file_header()
    }

    fn alignment_records<'a>(
        &'a mut self,
        header: &'a sam::Header,
    ) -> Box<dyn Iterator<Item = io::Result<sam::alignment::Record>> + 'a> {
        Box::new(
            self.records(header)
                .map(|result| result.and_then(|record| record.try_into_alignment_record(header))),
        )
    }
}

fn read_magic_number<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    let mut buf = [0; 4];
    reader.read_exact(&mut buf)?;

    if buf == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid CRAM header",
        ))
    }
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

    #[test]
    fn test_read_file_definition() -> Result<(), Box<dyn std::error::Error>> {
        let data = [
            0x43, 0x52, 0x41, 0x4d, // magic number = b"CRAM"
            0x03, 0x00, // format version = (3, 0)
            0x00, 0x68, 0xac, 0xf3, 0x06, 0x4d, 0xaa, 0x1e, 0x29, 0xa4, 0xa0, 0x8c, 0x56, 0xee,
            0x91, 0x9b, 0x91, 0x04, 0x21, 0x1f, // file ID
        ];

        let mut reader = Reader::new(&data[..]);
        let actual = reader.read_file_definition()?;

        let file_id = <[u8; 20]>::try_from(&data[6..])?;
        let expected = FileDefinition::new(Version::new(3, 0), file_id);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_magic_number() {
        let data = b"CRAM";
        let mut reader = &data[..];
        assert!(read_magic_number(&mut reader).is_ok());
    }

    #[test]
    fn test_read_magic_number_with_invalid_input() {
        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic_number(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof,
        ));

        let data = b"BAM\x01";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic_number(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData,
        ));
    }
}
