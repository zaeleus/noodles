//! CRAM writer.

pub(crate) mod builder;
mod collections;
pub(crate) mod container;
pub(crate) mod header;
pub(crate) mod num;
mod options;
pub(crate) mod record;

use std::io::{self, Write};

use noodles_fasta as fasta;
use noodles_sam::{self as sam, header::ReferenceSequences};

pub use self::builder::Builder;
use self::{
    container::write_container,
    header::{write_file_definition, write_file_header},
};
pub(crate) use self::{options::Options, record::Record};
use crate::{calculate_normalized_sequence_digest, FileDefinition};

const DEFAULT_SLICES_PER_CONTAINER: usize = 1;
const DEFAULT_RECORDS_PER_SLICE: usize = 10240;
pub(crate) const RECORDS_PER_CONTAINER: usize =
    DEFAULT_SLICES_PER_CONTAINER * DEFAULT_RECORDS_PER_SLICE;

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
    reference_sequence_repository: fasta::Repository,
    options: Options,
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
        let file_definition = FileDefinition::new(self.options.version, Default::default());
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
        let mut header = header.clone();

        add_missing_reference_sequence_checksums(
            &self.reference_sequence_repository,
            header.reference_sequences_mut(),
        )?;

        write_file_header(&mut self.inner, &header)
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
        self.write_file_definition()?;
        self.write_file_header(header)
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
    /// let record = sam::Record::default();
    /// writer.write_alignment_record(&header, &record)?;
    ///
    /// writer.try_finish(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, header: &sam::Header, record: Record) -> io::Result<()> {
        self.records.push(record);

        if self.records.len() >= self.records.capacity() {
            self.flush(header)?;
        }

        Ok(())
    }

    fn flush(&mut self, header: &sam::Header) -> io::Result<()> {
        write_container(
            &mut self.inner,
            &self.reference_sequence_repository,
            &self.options,
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
        let r = Record::try_from_alignment_record(header, record)?;
        self.write_record(header, r)
    }

    fn finish(&mut self, header: &sam::Header) -> io::Result<()> {
        self.try_finish(header)
    }
}

pub(crate) fn add_missing_reference_sequence_checksums(
    reference_sequence_repository: &fasta::Repository,
    reference_sequences: &mut ReferenceSequences,
) -> io::Result<()> {
    use indexmap::map::Entry;
    use sam::header::record::value::map::reference_sequence::{tag, Md5Checksum};

    for (name, reference_sequence) in reference_sequences {
        if let Entry::Vacant(entry) = reference_sequence
            .other_fields_mut()
            .entry(tag::MD5_CHECKSUM)
        {
            let sequence = reference_sequence_repository
                .get(name)
                .transpose()?
                .expect("missing reference sequence");

            let checksum = calculate_normalized_sequence_digest(&sequence[..]);

            entry.insert(Md5Checksum::from(checksum).to_string().into());
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use bstr::BString;

    use super::*;

    #[test]
    fn test_add_missing_reference_sequence_checksums() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZeroUsize;

        use fasta::record::{Definition, Sequence};
        use sam::header::record::value::{
            map::{reference_sequence::tag, ReferenceSequence},
            Map,
        };

        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        const SQ1_LN: NonZeroUsize = match NonZeroUsize::new(13) {
            Some(length) => length,
            None => unreachable!(),
        };

        let reference_sequences = vec![
            fasta::Record::new(
                Definition::new("sq0", None),
                Sequence::from(b"TTCACCCA".to_vec()),
            ),
            fasta::Record::new(
                Definition::new("sq1", None),
                Sequence::from(b"GATCTTACTTTTT".to_vec()),
            ),
        ];

        let sq0_md5_checksum = BString::from("be19336b7e15968f7ac7dc82493d9cd8");
        let sq1_md5_checksum = BString::from("d80f22a19aeeb623b3e4f746c762f21d");

        let repository = fasta::Repository::new(reference_sequences);

        let mut header = sam::Header::builder()
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
            .add_reference_sequence(
                "sq1",
                Map::<ReferenceSequence>::builder()
                    .set_length(SQ1_LN)
                    .insert(tag::MD5_CHECKSUM, sq1_md5_checksum.clone())
                    .build()?,
            )
            .build();

        add_missing_reference_sequence_checksums(&repository, header.reference_sequences_mut())?;

        let sq0 = header.reference_sequences().get(&b"sq0"[..]);
        assert_eq!(
            sq0.and_then(|rs| rs.other_fields().get(&tag::MD5_CHECKSUM)),
            Some(&sq0_md5_checksum)
        );

        let sq1 = header.reference_sequences().get(&b"sq1"[..]);
        assert_eq!(
            sq1.and_then(|rs| rs.other_fields().get(&tag::MD5_CHECKSUM)),
            Some(&sq1_md5_checksum)
        );

        Ok(())
    }
}
