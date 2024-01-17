//! CRAM writer.

pub(crate) mod builder;
pub(crate) mod container;
pub(crate) mod data_container;
pub(crate) mod header_container;
pub(crate) mod num;
mod options;
pub(crate) mod record;

pub use self::builder::Builder;
pub(crate) use self::options::Options;

use std::{
    io::{self, Write},
    mem,
};

use noodles_fasta as fasta;
use noodles_sam::{self as sam, header::ReferenceSequences};

use crate::{file_definition::Version, DataContainer, FileDefinition, Record, MAGIC_NUMBER};

/// A CRAM writer.
///
/// A call to [`Self::try_finish`] must be made before the writer is dropped.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use noodles_cram as cram;
/// use noodles_sam as sam;
///
/// let mut writer = cram::io::Writer::new(Vec::new());
///
/// let header = sam::Header::default();
/// writer.write_header(&header)?;
///
/// let record = cram::Record::default();
/// writer.write_record(&header, record)?;
///
/// writer.try_finish(&header)?;
/// # Ok::<(), io::Error>(())
/// ```
#[derive(Debug)]
pub struct Writer<W>
where
    W: Write,
{
    inner: W,
    reference_sequence_repository: fasta::Repository,
    options: Options,
    data_container_builder: crate::data_container::Builder,
    record_counter: u64,
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
    /// use noodles_cram as cram;
    /// let writer = cram::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Builder::default().build_with_writer(inner)
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let writer = cram::io::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
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
    /// let mut writer = cram::io::Writer::new(Vec::new());
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
    ///
    /// let mut writer = cram::io::Writer::new(Vec::new());
    /// writer.write_file_definition()?;
    ///
    /// assert_eq!(writer.get_ref(), &[
    ///     // magic number (CRAM)
    ///     0x43, 0x52, 0x41, 0x4d,
    ///     // format (major, minor)
    ///     0x03, 0x00,
    ///     // file ID
    ///     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ///     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    /// ]);
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
    /// let mut writer = cram::io::Writer::new(Vec::new());
    /// writer.write_file_definition()?;
    ///
    /// let header = sam::Header::default();
    /// writer.write_file_header(&header)?;
    ///
    /// writer.try_finish(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_file_header(&mut self, header: &sam::Header) -> io::Result<()> {
        use self::header_container::write_header_container;

        let mut header = header.clone();

        add_missing_reference_sequence_checksums(
            &self.reference_sequence_repository,
            header.reference_sequences_mut(),
        )?;

        write_header_container(&mut self.inner, &header)
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
    /// let header = sam::Header::builder().add_comment("noodles-cram").build();
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
    /// use noodles_sam as sam;
    ///
    /// let mut writer = cram::io::Writer::new(Vec::new());
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header)?;
    ///
    /// let record = cram::Record::default();
    /// writer.write_record(&header, record)?;
    ///
    /// writer.try_finish(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, header: &sam::Header, mut record: Record) -> io::Result<()> {
        use crate::data_container::builder::AddRecordError;

        loop {
            match self.data_container_builder.add_record(record) {
                Ok(_) => {
                    self.record_counter += 1;
                    return Ok(());
                }
                Err(e) => match e {
                    AddRecordError::ContainerFull(r) => {
                        record = r;
                        self.flush(header)?;
                    }
                    AddRecordError::SliceFull(r) => {
                        record = r;
                    }
                    _ => return Err(io::Error::from(io::ErrorKind::InvalidInput)),
                },
            }
        }
    }

    fn flush(&mut self, header: &sam::Header) -> io::Result<()> {
        use self::data_container::write_data_container;

        if self.data_container_builder.is_empty() {
            return Ok(());
        }

        let data_container_builder = mem::replace(
            &mut self.data_container_builder,
            DataContainer::builder(self.record_counter),
        );

        let base_count = data_container_builder.base_count();

        let data_container = data_container_builder.build(
            &self.options,
            &self.reference_sequence_repository,
            header,
        )?;

        write_data_container(&mut self.inner, &data_container, base_count)
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

fn write_file_definition<W>(writer: &mut W, file_definition: &FileDefinition) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(MAGIC_NUMBER)?;
    write_format(writer, file_definition.version())?;
    writer.write_all(file_definition.file_id())?;
    Ok(())
}

fn write_format<W>(writer: &mut W, version: Version) -> io::Result<()>
where
    W: Write,
{
    let format = [version.major(), version.minor()];
    writer.write_all(&format)
}

pub(crate) fn add_missing_reference_sequence_checksums(
    reference_sequence_repository: &fasta::Repository,
    reference_sequences: &mut ReferenceSequences,
) -> io::Result<()> {
    use indexmap::map::Entry;
    use sam::header::record::value::map::reference_sequence::{tag, Md5Checksum};

    use crate::data_container::slice::builder::calculate_normalized_sequence_digest;

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
