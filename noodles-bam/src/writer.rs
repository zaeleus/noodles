//! BAM writer.

pub mod record;

use std::{
    ffi::CString,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_bgzf as bgzf;
use noodles_sam::{
    self as sam,
    alignment::Record,
    header::{
        record::value::{map::ReferenceSequence, Map},
        ReferenceSequences,
    },
};

use self::record::encode_record;

/// A BAM writer.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use noodles_bam as bam;
/// use noodles_sam::{self as sam, alignment::Record};
///
/// let mut writer = bam::Writer::new(Vec::new());
///
/// let header = sam::Header::default();
/// writer.write_header(&header)?;
/// writer.write_reference_sequences(header.reference_sequences())?;
///
/// let record = Record::default();
/// writer.write_record(&header, &record)?;
/// # Ok::<(), io::Error>(())
/// ```
pub struct Writer<W> {
    inner: W,
    buf: Vec<u8>,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let writer = bam::Writer::from(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let mut writer = bam::Writer::from(Vec::new());
    /// assert!(writer.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let writer = bam::Writer::from(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
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
    /// use noodles_sam::{
    ///     self as sam,
    ///     header::record::value::{map::ReferenceSequence, Map},
    /// };
    ///
    /// let mut writer = bam::Writer::new(Vec::new());
    ///
    /// let header = sam::Header::builder()
    ///     .add_reference_sequence(Map::<ReferenceSequence>::new("sq0".parse()?, 8)?)
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
    /// use noodles_sam::{self as sam, alignment::Record};
    ///
    /// let mut writer = bam::Writer::new(Vec::new());
    ///
    /// let header = sam::Header::default();
    /// let record = Record::default();
    /// writer.write_record(&header, &record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, header: &sam::Header, record: &Record) -> io::Result<()> {
        self.buf.clear();
        encode_record(&mut self.buf, header, record)?;

        let block_size = u32::try_from(self.buf.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        self.inner.write_u32::<LittleEndian>(block_size)?;

        self.inner.write_all(&self.buf)?;

        Ok(())
    }
}

impl<W> Writer<bgzf::Writer<W>>
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
        Self::from(bgzf::Writer::new(writer))
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
}

impl<W> From<W> for Writer<W> {
    fn from(inner: W) -> Self {
        Self {
            inner,
            buf: Vec::new(),
        }
    }
}

impl<W> sam::AlignmentWriter for Writer<W>
where
    W: Write,
{
    fn write_alignment_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.write_header(header)?;
        self.write_reference_sequences(header.reference_sequences())?;
        Ok(())
    }

    fn write_alignment_record(&mut self, header: &sam::Header, record: &Record) -> io::Result<()> {
        self.write_record(header, record)
    }

    fn finish(&mut self, _: &sam::Header) -> io::Result<()> {
        Ok(())
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
    reference_sequence: &Map<ReferenceSequence>,
) -> io::Result<()>
where
    W: Write,
{
    let c_name = CString::new(reference_sequence.name().as_bytes())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let name = c_name.as_bytes_with_nul();

    let l_name =
        u32::try_from(name.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(l_name)?;
    writer.write_all(name)?;

    let l_ref = i32::try_from(usize::from(reference_sequence.length()))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(l_ref)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_sam::{record::Data, AlignmentWriter};

    use super::*;
    use crate::Reader;

    #[test]
    fn test_write_header() -> Result<(), Box<dyn std::error::Error>> {
        use sam::header::record::value::{
            map::{self, header::Version},
            Map,
        };

        let header = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
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
        use sam::header::record::value::map::reference_sequence::Name;

        let reference_sequences = [("sq0".parse()?, 8)]
            .into_iter()
            .map(|(name, len): (Name, usize)| {
                let sn = name.to_string();
                Map::<ReferenceSequence>::new(name, len).map(|rs| (sn, rs))
            })
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
    fn test_write_alignment_record() -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let record = Record::default();
        writer.write_alignment_record(&header, &record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().get_ref().as_slice());

        let mut record = Record::default();
        reader.read_record(&mut record)?;

        assert!(record.read_name().is_none());
        assert_eq!(record.flags(), sam::record::Flags::UNMAPPED);
        assert!(record.reference_sequence_id().is_none());
        assert!(record.alignment_start().is_none());
        assert!(record.mapping_quality().is_none());
        assert!(record.cigar().is_empty());
        assert!(record.mate_reference_sequence_id().is_none());
        assert!(record.mate_alignment_start().is_none());
        assert_eq!(record.template_length(), 0);
        assert!(record.sequence().is_empty());
        assert!(record.quality_scores().is_empty());
        assert!(record.data().is_empty());

        Ok(())
    }

    #[test]
    fn test_write_alignment_record_with_sequence_length_less_than_quality_scores_length(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();

        let mut record = Record::builder().set_sequence("AT".parse()?).build();
        *record.quality_scores_mut() = "NDLS".parse()?;

        assert!(writer.write_alignment_record(&header, &record).is_err());

        Ok(())
    }

    #[test]
    fn test_write_alignment_record_with_sequence_length_greater_than_quality_scores_length(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();

        let mut record = Record::builder().set_sequence("ATCG".parse()?).build();
        *record.quality_scores_mut() = "ND".parse()?;

        assert!(writer.write_alignment_record(&header, &record).is_err());

        Ok(())
    }

    #[test]
    fn test_write_alignment_record_with_no_sequence_and_with_quality_scores(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let mut record = Record::default();
        *record.quality_scores_mut() = "NDLS".parse()?;

        assert!(writer.write_alignment_record(&header, &record).is_err());

        Ok(())
    }

    #[test]
    fn test_write_alignment_record_with_sequence_and_no_quality_scores(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let record = Record::builder().set_sequence("ATCG".parse()?).build();

        writer.write_alignment_record(&header, &record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().get_ref().as_slice());

        let mut record = Record::default();
        reader.read_record(&mut record)?;

        let expected = "ATCG".parse()?;
        assert_eq!(record.sequence(), &expected);

        assert!(record.quality_scores().is_empty());

        Ok(())
    }

    #[test]
    fn test_write_alignment_record_with_sequence_and_quality_scores(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let sam_record = Record::builder()
            .set_sequence("ATCG".parse()?)
            .set_quality_scores("NDLS".parse()?)
            .build();

        writer.write_alignment_record(&header, &sam_record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().get_ref().as_slice());

        let mut record = Record::default();
        reader.read_record(&mut record)?;

        let expected = "ATCG".parse()?;
        assert_eq!(record.sequence(), &expected);

        assert_eq!(record.quality_scores(), sam_record.quality_scores());

        Ok(())
    }

    #[test]
    fn test_write_alignment_record_with_data() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_sam::record::data::{
            field::{Tag, Value},
            Field,
        };

        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let sam_record = Record::builder()
            .set_data(Data::try_from(vec![
                Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
                Field::new(Tag::AlignmentHitCount, Value::UInt8(1)),
            ])?)
            .build();

        writer.write_alignment_record(&header, &sam_record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().get_ref().as_slice());

        let mut record = Record::default();
        reader.read_record(&mut record)?;

        let bam_data = record.data();
        let mut fields = bam_data.values();

        assert_eq!(
            fields.next(),
            Some(&Field::new(
                Tag::ReadGroup,
                Value::String(String::from("rg0"))
            ))
        );

        assert_eq!(
            fields.next(),
            Some(&Field::new(Tag::AlignmentHitCount, Value::UInt8(1)))
        );

        assert!(fields.next().is_none());

        Ok(())
    }

    #[test]
    fn test_write_reference_sequence() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();
        let reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 8)?;
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
