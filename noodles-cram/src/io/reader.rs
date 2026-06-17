//! CRAM reader and record iterator.

mod builder;
pub(crate) mod collections;
pub(crate) mod container;
pub mod header;
pub(crate) mod num;
mod query;
mod records;

use std::{
    io::{self, Read, Seek, SeekFrom},
    vec,
};

use noodles_core::Region;
use noodles_fasta as fasta;
use noodles_sam as sam;

pub use self::{builder::Builder, container::Container, query::Query, records::Records};
use self::{container::read_container, header::read_header, records::decode_container_records};
use crate::{FileDefinition, crai};

/// A CRAM reader.
///
/// The CRAM format is comprised of four main parts: 1) a file definition, 2) a file header, 3) a
/// list of containers, and 4) an end-of-file (EOF) container.
///
/// # Examples
///
/// ```no_run
/// # use std::{fs::File, io};
/// use noodles_cram as cram;
/// use noodles_fasta as fasta;
///
/// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
/// let header = reader.read_header()?;
///
/// for result in reader.records(&header) {
///     let record = result?;
///     // ...
/// }
///
/// # Ok::<_, io::Error>(())
/// ```
pub struct Reader<R> {
    inner: R,
    reference_sequence_repository: fasta::Repository,
    container: Container,
    records: vec::IntoIter<sam::alignment::RecordBuf>,
    is_eof: bool,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let reader = cram::io::Reader::new(io::empty());
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let mut reader = cram::io::Reader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Unwraps and returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let reader = cram::io::Reader::new(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
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
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(inner: R) -> Self {
        Builder::default().build_from_reader(inner)
    }

    pub(crate) fn reference_sequence_repository(&self) -> &fasta::Repository {
        &self.reference_sequence_repository
    }

    /// Returns a CRAM header reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io::Read};
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
    ///
    /// let _raw_header = {
    ///     let mut raw_sam_header_reader = container_reader.raw_sam_header_reader()?;
    ///     let mut raw_header = String::new();
    ///     raw_sam_header_reader.read_to_string(&mut raw_header)?;
    ///     raw_sam_header_reader.discard_to_end()?;
    ///     raw_header
    /// };
    ///
    /// container_reader.discard_to_end()?;
    /// Ok::<_, std::io::Error>(())
    /// ```
    pub fn header_reader(&mut self) -> header::Reader<&mut R> {
        header::Reader::new(&mut self.inner)
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
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// let file_definition = reader.read_file_definition()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_file_definition(&mut self) -> io::Result<FileDefinition> {
        header::read_file_definition(&mut self.inner)
    }

    /// Reads the SAM header.
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
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// reader.read_file_definition()?;
    ///
    /// let header = reader.read_file_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_file_header(&mut self) -> io::Result<sam::Header> {
        header::read_file_header(&mut self.inner)
    }

    /// Reads the SAM header.
    ///
    /// This verifies the CRAM magic number, discards the file definition, and reads and parses the
    /// file header as a SAM header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// let header = reader.read_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<sam::Header> {
        read_header(&mut self.inner)
    }

    /// Reads a container.
    ///
    /// This returns `None` if the container header is the EOF container header, which signals the
    /// end of the stream.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram::{self as cram, io::reader::Container};
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// reader.read_header()?;
    ///
    /// let mut container = Container::default();
    ///
    /// while reader.read_container(&mut container)? != 0 {
    ///     // ...
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_container(&mut self, container: &mut Container) -> io::Result<usize> {
        read_container(&mut self.inner, container)
    }

    /// Reads the next alignment record.
    ///
    /// The record is read into `record`, replacing its previous contents. This advances through
    /// containers and slices as needed, decoding a whole container at a time and buffering its
    /// records internally so that subsequent calls are served from the buffer.
    ///
    /// The stream is expected to be at the start of a container.
    ///
    /// This returns the number of records read, which is `0` if the end of the stream was reached.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    /// use noodles_sam::alignment::RecordBuf;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// let header = reader.read_header()?;
    ///
    /// let mut record = RecordBuf::default();
    ///
    /// while reader.read_record_buf(&header, &mut record)? != 0 {
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_record_buf(
        &mut self,
        header: &sam::Header,
        record: &mut sam::alignment::RecordBuf,
    ) -> io::Result<usize> {
        loop {
            if let Some(r) = self.records.next() {
                *record = r;
                return Ok(1);
            }

            if !self.read_next_container_records(header)? {
                return Ok(0);
            }
        }
    }

    fn read_next_container_records(&mut self, header: &sam::Header) -> io::Result<bool> {
        if self.is_eof {
            return Ok(false);
        }

        if read_container(&mut self.inner, &mut self.container)? == 0 {
            self.is_eof = true;
            return Ok(false);
        }

        self.records =
            decode_container_records(&self.reference_sequence_repository, header, &self.container)?
                .into_iter();

        Ok(true)
    }

    /// Returns a iterator over records starting from the current stream position.
    ///
    /// The stream is expected to be at the start of a container.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    /// use noodles_fasta as fasta;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// let header = reader.read_header()?;
    ///
    /// for result in reader.records(&header) {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records<'r, 'h: 'r>(&'r mut self, header: &'h sam::Header) -> Records<'r, 'h, R> {
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
    /// This discards any records buffered by [`Self::read_record_buf`], as the stream position no
    /// longer matches the buffer.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io::{self, SeekFrom};
    /// use noodles_cram as cram;
    /// let mut reader = cram::io::Reader::new(io::empty());
    /// reader.seek(SeekFrom::Start(0))?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        let position = self.inner.seek(pos)?;

        self.records = Vec::new().into_iter();
        self.is_eof = false;

        Ok(position)
    }

    /// Returns the current position of the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let mut reader = cram::io::Reader::new(io::empty());
    /// let position = reader.position()?;
    /// assert_eq!(position, 0);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn position(&mut self) -> io::Result<u64> {
        self.inner.stream_position()
    }

    /// Returns an iterator over records that intersects the given region.
    ///
    /// To query for unmapped records, use [`Self::query_unmapped`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram::{self as cram, crai};
    /// use noodles_fasta as fasta;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let header = reader.read_header()?;
    /// let index = crai::fs::read("sample.cram.crai")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<'r, 'h: 'r, 'i: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        index: &'i crai::Index,
        region: &Region,
    ) -> io::Result<Query<'r, 'h, 'i, R>> {
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

    /// Returns an iterator of unmapped records after querying for the unmapped region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram::{self as cram, crai};
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let header = reader.read_header()?;
    /// let index = crai::fs::read("sample.cram.crai")?;
    /// let query = reader.query_unmapped(&header, &index)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn query_unmapped<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        index: &crai::Index,
    ) -> io::Result<impl Iterator<Item = io::Result<sam::alignment::RecordBuf>> + use<'r, 'h, R>>
    {
        let offset = index
            .iter()
            .find(|record| record.reference_sequence_id().is_none())
            .map(|record| SeekFrom::Start(record.offset()))
            .unwrap_or(SeekFrom::End(0));

        self.seek(offset)?;

        Ok(self.records(header).filter_map(|result| match result {
            Ok(record) => {
                if record.flags().is_unmapped() {
                    Some(Ok(record))
                } else {
                    None
                }
            }
            Err(e) => Some(Err(e)),
        }))
    }
}

impl<R> sam::alignment::io::Read<R> for Reader<R>
where
    R: Read,
{
    fn read_alignment_header(&mut self) -> io::Result<sam::Header> {
        self.read_header()
    }

    fn alignment_records<'a>(
        &'a mut self,
        header: &'a sam::Header,
    ) -> Box<dyn Iterator<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'a> {
        Box::new(
            self.records(header).map(|result| {
                result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
            }),
        )
    }
}

#[cfg(test)]
mod tests {
    use std::{io::Cursor, num::NonZero};

    use noodles_core::Position;
    use noodles_sam::{
        alignment::{
            RecordBuf,
            io::Write,
            record::{
                Flags,
                cigar::{Op, op::Kind},
            },
            record_buf::Sequence,
        },
        header::record::value::{Map, map::ReferenceSequence},
    };

    use super::*;
    use crate::io::writer;

    const REFERENCE_SEQUENCE: &[u8] = b"TTCACCCA";

    fn mapped_record(alignment_start: usize, bases: &[u8]) -> RecordBuf {
        RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(alignment_start).unwrap())
            .set_cigar([Op::new(Kind::Match, bases.len())].into_iter().collect())
            .set_sequence(Sequence::from(bases.to_vec()))
            .set_quality_scores(vec![45; bases.len()].into_iter().collect())
            .build()
    }

    fn write_fixture(
        records: &[RecordBuf],
    ) -> io::Result<(sam::Header, fasta::Repository, Vec<u8>)> {
        let reference_sequences = vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            fasta::record::Sequence::from(REFERENCE_SEQUENCE.to_vec()),
        )];

        let length = NonZero::try_from(REFERENCE_SEQUENCE.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let header = sam::Header::builder()
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(length))
            .build();

        let repository = fasta::Repository::new(reference_sequences);

        let mut buf = Vec::new();
        let mut writer = writer::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_writer(&mut buf);

        writer.write_header(&header)?;

        for record in records {
            writer.write_alignment_record(&header, record)?;
        }

        writer.try_finish(&header)?;

        Ok((header, repository, buf))
    }

    #[test]
    fn test_read_record_buf() -> io::Result<()> {
        let (header, repository, buf) = write_fixture(&[mapped_record(1, b"TTCA")])?;

        let mut reader = Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(&buf[..]);
        reader.read_header()?;

        let mut record = RecordBuf::default();

        assert_eq!(reader.read_record_buf(&header, &mut record)?, 1);
        assert_eq!(record.sequence(), &Sequence::from(b"TTCA"));
        assert_eq!(record.alignment_start(), Position::new(1));

        // Once the end-of-file container is reached, further reads return 0 idempotently rather
        // than attempting to read past it.
        assert_eq!(reader.read_record_buf(&header, &mut record)?, 0);
        assert_eq!(reader.read_record_buf(&header, &mut record)?, 0);

        Ok(())
    }

    #[test]
    fn test_read_record_buf_across_containers() -> io::Result<()> {
        // The writer emits a container per 10240 records, so one more than that spans two
        // containers and exercises advancing across the boundary.
        const RECORD_COUNT: usize = 10240 + 1;

        let records = vec![mapped_record(1, b"TTCA"); RECORD_COUNT];
        let (header, repository, buf) = write_fixture(&records)?;

        let mut reader = Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(&buf[..]);
        reader.read_header()?;

        let mut record = RecordBuf::default();
        let mut n = 0;

        while reader.read_record_buf(&header, &mut record)? != 0 {
            n += 1;
        }

        assert_eq!(n, RECORD_COUNT);

        Ok(())
    }

    #[test]
    fn test_read_record_buf_after_seek_discards_buffer() -> io::Result<()> {
        let records = [mapped_record(1, b"TTCA"), mapped_record(2, b"TCAC")];
        let (header, repository, buf) = write_fixture(&records)?;

        let mut reader = Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(Cursor::new(buf));
        reader.read_header()?;

        // Both records decode from a single container, so reading the first buffers the second
        // inside the reader.
        let container_position = reader.position()?;

        let mut record = RecordBuf::default();
        assert_eq!(reader.read_record_buf(&header, &mut record)?, 1);
        assert_eq!(record.alignment_start(), Position::new(1));

        // Seeking back to the container discards that buffered record, so the next read starts
        // over from the first record instead of yielding the stale buffer.
        reader.seek(SeekFrom::Start(container_position))?;

        assert_eq!(reader.read_record_buf(&header, &mut record)?, 1);
        assert_eq!(record.alignment_start(), Position::new(1));

        Ok(())
    }
}
