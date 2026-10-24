//! CRAM writer.

pub(crate) mod builder;
mod collections;
pub(crate) mod container;
mod context;
pub(crate) mod header;
pub(crate) mod num;
pub(crate) mod record;

use std::io::{self, Write};

use noodles_sam::{self as sam, alignment::io::Write as _};

pub use self::builder::Builder;
use self::{
    container::write_container,
    header::{write_file_definition, write_file_header, write_header},
};
pub(crate) use self::{context::Context, record::Record};
use crate::FileDefinition;

pub(crate) const DEFAULT_SLICES_PER_CONTAINER: usize = 1;

/// A CRAM writer.
///
/// A call to [`Self::try_finish`] must be made before the writer is dropped.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use noodles_cram as cram;
/// use noodles_sam::{self as sam, alignment::io::Write};
///
/// let mut writer = cram::io::Writer::new(io::sink());
///
/// let header = sam::Header::default();
/// writer.write_header(&header)?;
///
/// let record = sam::Record::default();
/// writer.write_alignment_record(&header, &record)?;
///
/// writer.try_finish(&header)?;
/// # Ok::<(), io::Error>(())
/// ```
#[derive(Debug)]
pub struct Writer<W> {
    inner: W,
    context: Context,
    records: Vec<Record>,
    record_counter: u64,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let writer = cram::io::Writer::new(io::sink());
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
    /// use noodles_cram as cram;
    /// let mut writer = cram::io::Writer::new(io::sink());
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
    /// use noodles_cram as cram;
    /// let mut writer = cram::io::Writer::new(io::sink());
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
    /// Creates a new CRAM writer with default options.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let writer = cram::io::Writer::new(io::sink());
    /// ```
    pub fn new(inner: W) -> Self {
        Builder::default().build_from_writer(inner)
    }

    /// Attempts to finish the output stream by writing any pending containers and a final EOF
    /// container.
    ///
    /// This is typically only manually called if the underlying stream is needed before the writer
    /// is dropped.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// use noodles_sam as sam;
    ///
    /// let header = sam::Header::default();
    /// let mut writer = cram::io::Writer::new(io::sink());
    /// writer.try_finish(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_finish(&mut self, header: &sam::Header) -> io::Result<()> {
        use self::container::write_eof_container;
        self.flush(header)?;
        write_eof_container(&mut self.inner)
    }

    /// Writes a CRAM file definition.
    ///
    /// The file ID is set as a blank value (`[0x00; 20]`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let mut writer = cram::io::Writer::new(io::sink());
    /// writer.write_file_definition()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_file_definition(&mut self) -> io::Result<()> {
        let file_definition = FileDefinition::new(self.context.version, Default::default());
        write_file_definition(&mut self.inner, &file_definition)
    }

    /// Writes a CRAM file header container.
    ///
    /// The position of the stream is expected to be directly after the file definition.
    ///
    /// Entries in the reference sequence dictionary that are missing MD5 checksums (`M5`) will
    /// automatically be calculated and added to the written record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// use noodles_sam as sam;
    ///
    /// let mut writer = cram::io::Writer::new(io::sink());
    /// writer.write_file_definition()?;
    ///
    /// let header = sam::Header::default();
    /// writer.write_file_header(&header)?;
    ///
    /// writer.try_finish(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_file_header(&mut self, header: &sam::Header) -> io::Result<()> {
        write_file_header(
            &mut self.inner,
            &self.context.reference_sequence_repository,
            header,
        )
    }

    /// Writes a SAM header.
    ///
    /// This writes the CRAM magic number, the file definition, and file header using the given SAM
    /// header.
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// use noodles_sam as sam;
    ///
    /// let mut writer = cram::io::Writer::new(io::sink());
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        let file_definition = FileDefinition::new(self.context.version, Default::default());

        write_header(
            &mut self.inner,
            &self.context.reference_sequence_repository,
            &file_definition,
            header,
        )
    }

    /// Writes a CRAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// use noodles_sam::{self as sam, alignment::io::Write};
    ///
    /// let mut writer = cram::io::Writer::new(io::sink());
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header)?;
    ///
    /// let record = cram::Record::default();
    /// writer.write_record(&header, &record)?;
    ///
    /// writer.try_finish(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(
        &mut self,
        header: &sam::Header,
        record: &crate::Record<'_>,
    ) -> io::Result<()> {
        self.write_alignment_record(header, record)
    }

    fn add_record(&mut self, header: &sam::Header, record: Record) -> io::Result<()> {
        self.records.push(record);

        if self.records.len() >= self.records.capacity() {
            self.flush(header)?;
        }

        Ok(())
    }

    fn flush(&mut self, header: &sam::Header) -> io::Result<()> {
        write_container(
            &mut self.inner,
            &self.context,
            header,
            self.record_counter,
            &mut self.records,
        )?;

        let record_count = u64::try_from(self.records.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        self.record_counter += record_count;

        self.records.clear();

        Ok(())
    }
}

impl<W> sam::alignment::io::Write for Writer<W>
where
    W: Write,
{
    fn write_alignment_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.write_header(header)
    }

    fn write_alignment_record(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
    ) -> io::Result<()> {
        let record = Record::try_from_alignment_record(
            &self.context.reference_sequence_repository,
            header,
            record,
        )?;

        self.add_record(header, record)
    }

    fn finish(&mut self, header: &sam::Header) -> io::Result<()> {
        self.try_finish(header)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use noodles_core::Position;
    use noodles_sam::alignment::{
        io::Write as _,
        record::{
            Flags,
            cigar::{Op, op::Kind},
        },
    };

    use super::*;

    fn multi_slice_fixture() -> io::Result<(Vec<u8>, fasta::Repository, crate::crai::Index)> {
        let repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            fasta::record::Sequence::from(vec![b'A'; 64]),
        )]);
        let header: sam::Header = "@SQ\tSN:sq0\tLN:64\n".parse().unwrap();
        let mut writer = Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_writer(Vec::new());
        writer.context.records_per_slice = 2;
        writer.write_header(&header)?;
        for i in 0..4 {
            let record = sam::alignment::RecordBuf::builder()
                .set_name(format!("r{i}"))
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::new(i + 1).unwrap())
                .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
                .set_sequence(b"AAAA".to_vec().into())
                .build();
            writer.write_alignment_record(&header, &record)?;
        }
        writer.try_finish(&header)?;
        let src = writer.into_inner();
        let mut reader = crate::io::Reader::new(Cursor::new(&src));
        reader.read_header()?;
        let offset = reader.get_ref().position();
        let mut container = crate::io::reader::Container::default();
        reader.read_container(&mut container)?;
        assert_eq!(container.header().landmarks().len(), 2);
        let mut index = Vec::new();
        for (i, (&landmark, slice)) in container
            .header()
            .landmarks()
            .iter()
            .zip(container.slices())
            .enumerate()
        {
            let slice = slice?;
            let crate::container::ReferenceSequenceContext::Some(context) =
                slice.header().reference_sequence_context()
            else {
                unreachable!()
            };
            index.push(crate::crai::Record::new(
                Some(context.reference_sequence_id()),
                Some(context.alignment_start()),
                context.alignment_span(),
                offset,
                landmark as u64,
                (container
                    .header()
                    .landmarks()
                    .get(i + 1)
                    .copied()
                    .unwrap_or(container.src.len())
                    - landmark) as u64,
            ));
        }
        Ok((src, repository, index))
    }

    #[test]
    fn test_query_does_not_repeat_records_from_other_slices() -> io::Result<()> {
        let (src, repository, index) = multi_slice_fixture()?;
        let mut reader = crate::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(Cursor::new(src));
        let header = reader.read_header()?;
        let region = "sq0:1-64".parse().unwrap();
        let records = reader
            .query(&header, &index, &region)?
            .records()
            .collect::<io::Result<Vec<_>>>()?;
        let names: Vec<_> = records
            .iter()
            .map(|r| r.name().unwrap().to_string())
            .collect();
        assert_eq!(names, ["r0", "r1", "r2", "r3"]);
        Ok(())
    }

    #[cfg(feature = "async")]
    #[tokio::test]
    async fn test_async_query_does_not_repeat_records_from_other_slices() -> io::Result<()> {
        use futures::TryStreamExt;
        let (src, repository, index) = multi_slice_fixture()?;
        let mut reader = crate::r#async::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(Cursor::new(src));
        let header = reader.read_header().await?;
        let region = "sq0:1-64".parse().unwrap();
        let records: Vec<_> = reader
            .query(&header, &index, &region)?
            .records()
            .try_collect()
            .await?;
        let names: Vec<_> = records
            .iter()
            .map(|r| r.name().unwrap().to_string())
            .collect();
        assert_eq!(names, ["r0", "r1", "r2", "r3"]);
        Ok(())
    }
}
