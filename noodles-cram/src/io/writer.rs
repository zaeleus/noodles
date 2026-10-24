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
mod reference_range_tests {
    use std::{
        io::{self, Cursor},
        num::NonZero,
    };

    use crate as cram;
    use noodles_core::Position;
    use noodles_fasta as fasta;
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::{
        self as sam,
        alignment::{
            RecordBuf,
            record::{
                Flags, MappingQuality,
                cigar::{Op, op::Kind},
            },
        },
    };

    type Sequences = [(&'static str, Vec<u8>)];

    fn header(sequences: &Sequences) -> sam::Header {
        let mut text = String::new();
        for (name, sequence) in sequences {
            text += &format!("@SQ\tSN:{name}\tLN:{}\n", sequence.len());
        }
        text.parse().unwrap()
    }

    fn repository(sequences: &Sequences) -> fasta::Repository {
        let record = |(name, sequence): &(&str, Vec<u8>)| {
            let definition = fasta::record::Definition::new(*name, None);
            fasta::Record::new(definition, sequence.clone().into())
        };
        fasta::Repository::new(sequences.iter().map(record).collect::<Vec<_>>())
    }

    /// An indexed FASTA repository, which reads ranges of sequences.
    fn indexed_repository(sequences: &Sequences) -> fasta::Repository {
        let (mut data, mut index) = (Vec::new(), Vec::new());
        for (name, sequence) in sequences {
            data.extend(format!(">{name}\n").bytes());
            let (offset, len) = (data.len() as u64, sequence.len() as u64);
            for line in sequence.chunks(60) {
                data.extend(line.iter().chain(b"\n"));
            }
            let (bases, width) = (NonZero::new(60).unwrap(), NonZero::new(61).unwrap());
            index.push(fasta::fai::Record::new(*name, len, offset, bases, width));
        }
        let reader = fasta::io::IndexedReader::new(Cursor::new(data), index.into());
        fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(reader))
    }

    /// Writes each group as its own container, without adding a writer API.
    fn write(sequences: &Sequences, groups: &[Vec<RecordBuf>]) -> io::Result<Vec<u8>> {
        let header = header(sequences);
        let mut writer = cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repository(sequences))
            .build_from_writer(Vec::new());
        writer.write_header(&header)?;
        for group in groups {
            for record in group {
                writer.write_alignment_record(&header, record)?;
            }
            writer.flush(&header)?;
        }
        writer.try_finish(&header)?;
        Ok(writer.into_inner())
    }

    fn read(src: &[u8], repository: fasta::Repository) -> io::Result<Vec<RecordBuf>> {
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(src);
        let header = reader.read_header()?;
        reader.records(&header).collect()
    }

    /// A mapped record whose matches copy the reference, except for a mismatch at each match start.
    fn record(name: &str, id: usize, start: usize, ops: &[(Kind, usize)], r: &[u8]) -> RecordBuf {
        let (mut sequence, mut position) = (Vec::new(), start - 1);
        for &(kind, n) in ops {
            match kind {
                Kind::Match => {
                    sequence.push(b'N');
                    sequence.extend(&r[position + 1..position + n]);
                    position += n;
                }
                Kind::Deletion | Kind::Skip => position += n,
                _ => sequence.extend(b"GGGG".iter().take(n)),
            }
        }
        let cigar = ops.iter().map(|&(kind, n)| Op::new(kind, n)).collect();
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::empty())
            .set_reference_sequence_id(id)
            .set_alignment_start(Position::new(start).unwrap())
            .set_mapping_quality(MappingQuality::new(30).unwrap())
            .set_cigar(cigar)
            .set_quality_scores(vec![30; sequence.len()].into())
            .set_sequence(sequence.into())
            .build()
    }

    #[test]
    fn indexed_reference_spans_match_whole_reference() -> io::Result<()> {
        use Kind::{Deletion, Insertion, Match, Skip, SoftClip};
        // `sq0` is long, so slices load only their spans. `sq1` is short, so it is loaded whole.
        let base = |i: usize| b"ACGTNacgtn"[i * 7 % 11 % 10];
        let sequences = [
            ("sq0", (0..20_000).map(base).collect::<Vec<_>>()),
            ("sq1", (0..200).map(base).collect()),
        ];
        let (sq0, sq1) = (&sequences[0].1, &sequences[1].1);
        #[rustfmt::skip]
        let complex = [(SoftClip, 2), (Match, 10), (Insertion, 1), (Match, 5), (Deletion, 2), (Match, 8)];
        let groups = [
            // A single-reference slice, which has a reference checksum.
            vec![
                record("a0", 0, 5000, &complex, sq0),
                record("a1", 0, 5010, &[(Match, 10), (Skip, 3), (Match, 10)], sq0),
            ],
            // A multi-reference slice.
            vec![
                record("b0", 0, 6000, &complex, sq0),
                record("b1", 1, 10, &[(Match, 30)], sq1),
                record("b2", 1, 150, &complex, sq1),
            ],
        ];
        let src = write(&sequences, &groups)?;
        let indexed = indexed_repository(&sequences);
        let actual = read(&src, indexed.clone())?;
        assert_eq!(actual, read(&src, repository(&sequences))?);
        // Only `sq1` was loaded whole.
        assert_eq!(indexed.len(), 1);
        indexed.get(b"sq1");
        assert_eq!(indexed.len(), 1);

        // A reference base that differs within a slice span fails the slice checksum either way.
        let mut corrupted = sequences.clone();
        corrupted[0].1[5004] = b'T';
        for repository in [repository(&corrupted), indexed_repository(&corrupted)] {
            let error = read(&src, repository).unwrap_err();
            assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        }
        Ok(())
    }
}
