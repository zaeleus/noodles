pub(crate) mod sam_record;

use std::{
    ffi::CString,
    io::{self, Write},
    mem,
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_bgzf as bgzf;
use noodles_sam::{
    self as sam,
    header::{ReferenceSequence, ReferenceSequences},
};

use super::Record;

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
    /// Creates a BAM writer with a default compression level.
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
        write_header(&mut self.inner, header)
    }

    /// Writes SAM reference sequences.
    ///
    /// The reference sequences here are typically the same as the reference sequences in the SAM
    /// header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{self as sam, header::ReferenceSequence};
    ///
    /// let mut writer = bam::Writer::new(Vec::new());
    ///
    /// let header = sam::Header::builder()
    ///     .add_reference_sequence(ReferenceSequence::new("sq0", 8)?)
    ///     .add_comment("noodles-bam")
    ///     .build();
    ///
    /// writer.write_header(&header)?;
    /// writer.write_reference_sequences(header.reference_sequences())?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write_reference_sequences(
        &mut self,
        reference_sequences: &ReferenceSequences,
    ) -> io::Result<()> {
        write_reference_sequences(&mut self.inner, reference_sequences)
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
        write_record(&mut self.inner, record)
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
        sam_record::write_sam_record(&mut self.inner, reference_sequences, record)
    }
}

fn write_header<W>(writer: &mut W, header: &sam::Header) -> io::Result<()>
where
    W: Write,
{
    use super::MAGIC_NUMBER;

    writer.write_all(MAGIC_NUMBER)?;

    let text = header.to_string();
    let l_text =
        i32::try_from(text.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(l_text)?;

    writer.write_all(text.as_bytes())?;

    Ok(())
}

fn write_reference_sequences<W>(
    writer: &mut W,
    reference_sequences: &ReferenceSequences,
) -> io::Result<()>
where
    W: Write,
{
    let n_ref = i32::try_from(reference_sequences.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(n_ref)?;

    for reference_sequence in reference_sequences.values() {
        write_reference_sequence(writer, reference_sequence)?;
    }

    Ok(())
}

fn write_reference_sequence<W>(
    writer: &mut W,
    reference_sequence: &ReferenceSequence,
) -> io::Result<()>
where
    W: Write,
{
    let c_name = CString::new(reference_sequence.name())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let name = c_name.as_bytes_with_nul();

    let l_name =
        u32::try_from(name.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(l_name)?;
    writer.write_all(name)?;

    let l_ref = reference_sequence.len();
    writer.write_i32::<LittleEndian>(l_ref)?;

    Ok(())
}

fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: Write,
{
    let block_size = u32::try_from(record.block_size())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(block_size)?;

    writer.write_i32::<LittleEndian>(record.ref_id)?;
    writer.write_i32::<LittleEndian>(record.pos)?;

    let l_read_name = u8::try_from(record.read_name.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u8(l_read_name)?;

    let mapq = u8::from(record.mapping_quality());
    writer.write_u8(mapq)?;

    writer.write_u16::<LittleEndian>(record.bin())?;

    let n_cigar_op = u16::try_from(record.cigar().len() / mem::size_of::<u32>())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u16::<LittleEndian>(n_cigar_op)?;

    let flag = u16::from(record.flags());
    writer.write_u16::<LittleEndian>(flag)?;

    let l_seq = u32::try_from(record.sequence().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(l_seq)?;

    writer.write_i32::<LittleEndian>(record.next_ref_id)?;
    writer.write_i32::<LittleEndian>(record.next_pos)?;

    writer.write_i32::<LittleEndian>(record.template_length())?;

    writer.write_all(&record.read_name)?;

    for &raw_op in record.cigar().iter() {
        writer.write_u32::<LittleEndian>(raw_op)?;
    }

    writer.write_all(record.sequence().as_ref())?;
    writer.write_all(record.quality_scores())?;
    writer.write_all(record.data().as_ref())?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_sam::record::Data;

    use crate::{record::sequence::Base, Reader, Record};

    use super::*;

    #[test]
    fn test_write_header() -> io::Result<()> {
        use sam::header::header::{Header, Version};

        let header = sam::Header::builder()
            .set_header(Header::new(Version::new(1, 6)))
            .build();

        let mut buf = Vec::new();
        write_header(&mut buf, &header)?;

        let mut expected = vec![
            b'B', b'A', b'M', 0x01, // magic
            0x0b, 0x00, 0x00, 0x00, // l_text = 11
        ];
        expected.extend_from_slice(b"@HD\tVN:1.6\n"); // text

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_reference_sequences() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequences = [("sq0", 8)]
            .into_iter()
            .map(|(name, len)| ReferenceSequence::new(name, len).map(|rs| (name.into(), rs)))
            .collect::<Result<_, _>>()?;

        let mut buf = Vec::new();
        write_reference_sequences(&mut buf, &reference_sequences)?;

        let expected = [
            0x01, 0x00, 0x00, 0x00, // n_ref = 1
            0x04, 0x00, 0x00, 0x00, // ref[0].l_name = 4
            b's', b'q', b'0', 0x00, // ref[0].name = b"sq0\x00"
            0x08, 0x00, 0x00, 0x00, // ref[0].l_ref = 8
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_record() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();
        let record = Record::default();
        write_record(&mut buf, &record)?;

        let expected = [
            0x22, 0x00, 0x00, 0x00, // block_size = 34
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x02, // l_read_name = 2
            0xff, // mapq = 255
            0x48, 0x12, // bin = 4680
            0x00, 0x00, // n_cigar_op = 0
            0x04, 0x00, // flag = 4
            0x00, 0x00, 0x00, 0x00, // l_seq = 0
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            0x2a, 0x00, // read_name = "*\x00"
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_sam_record() -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let sam_record = sam::Record::default();
        writer.write_sam_record(header.reference_sequences(), &sam_record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().as_slice());

        let mut record = Record::default();
        reader.read_record(&mut record)?;

        assert_eq!(record.read_name()?.to_bytes(), b"*");
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
        assert_eq!(record.template_length(), 0);
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

        let mut record = sam::Record::builder().set_sequence("AT".parse()?).build()?;
        *record.quality_scores_mut() = "NDLS".parse()?;

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

        let mut record = sam::Record::builder()
            .set_sequence("ATCG".parse()?)
            .build()?;
        *record.quality_scores_mut() = "ND".parse()?;

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
        let sam_record = sam::Record::builder()
            .set_sequence("ATCG".parse()?)
            .build()?;

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
        assert_eq!(**actual, expected);

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
            .build()?;

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
        assert_eq!(**actual, expected);

        Ok(())
    }

    #[test]
    fn test_write_sam_record_with_data() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_sam::record::data::{
            field::{Tag as SamTag, Value as SamValue},
            Field as SamField,
        };

        use crate::record::data::{field::Value, Field};

        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let sam_record = sam::Record::builder()
            .set_data(Data::try_from(vec![
                SamField::new(SamTag::ReadGroup, SamValue::String(String::from("rg0"))),
                SamField::new(SamTag::AlignmentHitCount, SamValue::Int(1)),
            ])?)
            .build()?;

        writer.write_sam_record(header.reference_sequences(), &sam_record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().as_slice());

        let mut record = Record::default();
        reader.read_record(&mut record)?;

        let bam_data = record.data();
        let mut fields = bam_data.values();

        assert_eq!(
            fields.next().transpose()?,
            Some(Field::new(
                SamTag::ReadGroup,
                Value::String(String::from("rg0"))
            ),)
        );

        assert_eq!(
            fields.next().transpose()?,
            Some(Field::new(SamTag::AlignmentHitCount, Value::UInt8(1)))
        );

        assert!(fields.next().is_none());

        Ok(())
    }

    #[test]
    fn test_write_reference_sequence() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();
        let reference_sequence = ReferenceSequence::new("sq0", 8)?;
        write_reference_sequence(&mut buf, &reference_sequence)?;

        let expected = [
            0x04, 0x00, 0x00, 0x00, // l_name = 4
            0x73, 0x71, 0x30, 0x00, // name = b"sq0\x00"
            0x08, 0x00, 0x00, 0x00, // l_ref = 8
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
