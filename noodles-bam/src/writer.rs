pub(crate) mod record;

use std::{
    ffi::CString,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_bgzf as bgzf;
use noodles_sam::{
    self as sam,
    header::{ReferenceSequence, ReferenceSequences},
};

use super::{Record, MAGIC_NUMBER};

/// A BAM writer.
///
/// Since the raw text header and `bam::Record` are immutable, BAM files are created by encoding a
/// SAM header and SAM records.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use noodles_bam as bam;
/// use noodles_sam as sam;
///
/// let mut writer = bam::Writer::new(Vec::new());
///
/// let header = sam::Header::builder().add_comment("noodles-bam").build();
/// writer.write_header(&header)?;
/// writer.write_reference_sequences(header.reference_sequences())?;
///
/// let record = sam::Record::default();
/// writer.write_sam_record(header.reference_sequences(), &record)?;
/// # Ok::<(), io::Error>(())
/// ```
pub struct Writer<W>
where
    W: Write,
{
    inner: bgzf::Writer<W>,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a new writer with a default compression level.
    ///
    /// The given stream is wrapped in a BGZF encoder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let writer = bam::Writer::new(Vec::new());
    /// ```
    pub fn new(writer: W) -> Self {
        Self {
            inner: bgzf::Writer::new(writer),
        }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let writer = bam::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        self.inner.get_ref()
    }

    /// Attempts to finish the output stream.
    ///
    /// This is typically only manually called if the underlying stream is needed before the writer
    /// is dropped.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    /// let mut writer = bam::Writer::new(Vec::new());
    /// writer.try_finish()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_finish(&mut self) -> io::Result<()> {
        self.inner.try_finish()
    }

    /// Writes a SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    ///
    /// let mut writer = bam::Writer::new(Vec::new());
    ///
    /// let header = sam::Header::builder().add_comment("noodles-bam").build();
    /// writer.write_header(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.inner.write_all(MAGIC_NUMBER)?;

        let text = header.to_string();
        let l_text = text.len() as i32;
        self.inner.write_i32::<LittleEndian>(l_text)?;

        self.inner.write_all(text.as_bytes())?;

        Ok(())
    }

    /// Writes SAM reference sequences.
    ///
    /// The reference sequences here are typically the same as the reference sequences in the SAM
    /// header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    ///
    /// let mut writer = bam::Writer::new(Vec::new());
    ///
    /// let header = sam::Header::builder()
    ///     .add_reference_sequence(sam::header::ReferenceSequence::new(String::from("sq0"), 8))
    ///     .add_comment("noodles-bam")
    ///     .build();
    ///
    /// writer.write_header(&header)?;
    /// writer.write_reference_sequences(header.reference_sequences())?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_reference_sequences(
        &mut self,
        reference_sequences: &ReferenceSequences,
    ) -> io::Result<()> {
        let n_ref = reference_sequences.len() as i32;
        self.inner.write_i32::<LittleEndian>(n_ref)?;

        for reference_sequence in reference_sequences.values() {
            write_reference(&mut self.inner, reference_sequence)?;
        }

        Ok(())
    }

    /// Writes a BAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    /// let mut writer = bam::Writer::new(Vec::new());
    /// let record = bam::Record::default();
    /// writer.write_record(&record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        let block_size = record.len() as u32;
        self.inner.write_u32::<LittleEndian>(block_size)?;
        self.inner.write_all(record)
    }

    /// Writes a SAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    ///
    /// let mut writer = bam::Writer::new(Vec::new());
    ///
    /// let reference_sequences = sam::header::ReferenceSequences::new();
    /// let record = sam::Record::default();
    /// writer.write_sam_record(&reference_sequences, &record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_sam_record(
        &mut self,
        reference_sequences: &ReferenceSequences,
        record: &sam::Record,
    ) -> io::Result<()> {
        record::write_sam_record(&mut self.inner, reference_sequences, record)
    }
}

impl<W> Drop for Writer<W>
where
    W: Write,
{
    fn drop(&mut self) {
        let _ = self.try_finish();
    }
}

fn write_reference<W>(writer: &mut W, reference_sequence: &ReferenceSequence) -> io::Result<()>
where
    W: Write,
{
    let c_name = CString::new(reference_sequence.name())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let name = c_name.as_bytes_with_nul();

    let l_name = name.len() as i32;
    writer.write_i32::<LittleEndian>(l_name)?;
    writer.write_all(name)?;

    let l_ref = reference_sequence.len() as i32;
    writer.write_i32::<LittleEndian>(l_ref)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_sam::record::Data;

    use crate::{record::sequence::Base, Reader, Record};

    use super::*;

    #[test]
    fn test_write_header() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::builder()
            .set_header(sam::header::header::Header::default())
            .build();

        writer.write_header(&header)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().as_slice());
        let actual = reader.read_header()?;

        let expected = "@HD\tVN:1.6\n";

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_write_reference_sequences() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::builder()
            .add_reference_sequence(sam::header::ReferenceSequence::new(String::from("sq0"), 8))
            .set_header(sam::header::header::Header::default())
            .build();

        writer.write_header(&header)?;
        writer.write_reference_sequences(header.reference_sequences())?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().as_slice());
        reader.read_header()?;
        let actual = reader.read_reference_sequences()?;

        assert_eq!(actual.len(), 1);
        assert_eq!(
            &actual[0],
            &sam::header::ReferenceSequence::new(String::from("sq0"), 8)
        );

        Ok(())
    }

    #[test]
    fn test_write_sam_record() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let sam_record = sam::Record::default();
        writer.write_sam_record(header.reference_sequences(), &sam_record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().as_slice());

        let mut record = Record::default();
        reader.read_record(&mut record)?;

        assert_eq!(record.read_name(), b"*\0");
        assert_eq!(record.flags(), sam::record::Flags::UNMAPPED);
        assert!(record.reference_sequence_id().is_none());
        assert!(record.position().is_none());
        assert_eq!(
            record.mapping_quality(),
            sam::record::MappingQuality::from(255)
        );
        assert!(record.cigar().is_empty());
        assert!(record.mate_reference_sequence_id().is_none());
        assert!(record.mate_position().is_none());
        assert_eq!(record.template_len(), 0);
        assert!(record.sequence().is_empty());
        assert!(record.quality_scores().is_empty());
        assert!(record.data().is_empty());

        Ok(())
    }

    #[test]
    fn test_write_sam_record_with_sequence_length_less_than_quality_scores_length(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let record = sam::Record::builder()
            .set_sequence("AT".parse()?)
            .set_quality_scores("NDLS".parse()?)
            .build();

        assert!(writer
            .write_sam_record(header.reference_sequences(), &record)
            .is_err());

        Ok(())
    }

    #[test]
    fn test_write_sam_record_with_sequence_length_greater_than_quality_scores_length(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let record = sam::Record::builder()
            .set_sequence("ATCG".parse()?)
            .set_quality_scores("ND".parse()?)
            .build();

        assert!(writer
            .write_sam_record(header.reference_sequences(), &record)
            .is_err());

        Ok(())
    }

    #[test]
    fn test_write_sam_record_with_sequence_and_no_quality_scores(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let sam_record = sam::Record::builder().set_sequence("ATCG".parse()?).build();
        writer.write_sam_record(header.reference_sequences(), &sam_record)?;

        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().as_slice());

        let mut record = Record::default();
        reader.read_record(&mut record)?;

        let actual: Vec<_> = record.sequence().bases().collect();
        let expected = [Base::A, Base::T, Base::C, Base::G];
        assert_eq!(actual, expected);

        let actual = record.quality_scores();
        let expected = [255, 255, 255, 255];
        assert_eq!(*actual, expected);

        Ok(())
    }

    #[test]
    fn test_write_sam_record_with_sequence_and_quality_scores(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let sam_record = sam::Record::builder()
            .set_sequence("ATCG".parse()?)
            .set_quality_scores("NDLS".parse()?)
            .build();

        writer.write_sam_record(header.reference_sequences(), &sam_record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().as_slice());

        let mut record = Record::default();
        reader.read_record(&mut record)?;

        let actual: Vec<_> = record.sequence().bases().collect();
        let expected = [Base::A, Base::T, Base::C, Base::G];
        assert_eq!(actual, expected);

        let actual = record.quality_scores();
        let expected = [45, 35, 43, 50];
        assert_eq!(*actual, expected);

        Ok(())
    }

    #[test]
    fn test_write_sam_record_with_data() -> io::Result<()> {
        use noodles_sam::record::data::{
            field::{Tag as SamTag, Value as SamValue},
            Field as SamField,
        };

        use crate::record::data::{field::Value, Field};

        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let sam_record = sam::Record::builder()
            .set_data(Data::from(vec![
                SamField::new(SamTag::ReadGroup, SamValue::String(String::from("rg0"))),
                SamField::new(SamTag::AlignmentHitCount, SamValue::Int32(1)),
            ]))
            .build();

        writer.write_sam_record(header.reference_sequences(), &sam_record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().as_slice());

        let mut record = Record::default();
        reader.read_record(&mut record)?;

        let bam_data = record.data();
        let mut fields = bam_data.fields();

        assert_eq!(
            fields.next().transpose()?,
            Some(Field::new(
                SamTag::ReadGroup,
                Value::String(String::from("rg0"))
            ),)
        );

        assert_eq!(
            fields.next().transpose()?,
            Some(Field::new(SamTag::AlignmentHitCount, Value::Int32(1)))
        );

        assert!(fields.next().is_none());

        Ok(())
    }
}
