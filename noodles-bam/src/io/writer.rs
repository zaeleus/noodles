//! BAM writer.

mod builder;
mod header;
pub(crate) mod num;

use std::io::{self, Write};

use noodles_bgzf as bgzf;
use noodles_sam::{self as sam, alignment::io::Write as _};

pub use self::builder::Builder;
use self::num::write_u32_le;
use crate::Record;

/// A BAM writer.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use noodles_bam as bam;
/// use noodles_sam as sam;
///
/// let mut writer = bam::io::Writer::new(io::sink());
///
/// let header = sam::Header::default();
/// writer.write_header(&header)?;
///
/// let record = bam::Record::default();
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
    /// let writer = bam::io::Writer::from(Vec::new());
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
    /// let mut writer = bam::io::Writer::from(Vec::new());
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
    /// let writer = bam::io::Writer::from(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }

    /// Writes a SAM header.
    ///
    /// This writes the BAM magic number, the raw SAM header, and a copy of the reference sequence
    /// dictionary as binary reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    ///
    /// let mut writer = bam::io::Writer::new(io::sink());
    ///
    /// let header = sam::Header::builder().add_comment("noodles-bam").build();
    /// writer.write_header(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        use self::header::write_header;
        write_header(&mut self.inner, header)
    }

    /// Writes a BAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    ///
    /// let header = sam::Header::default();
    ///
    /// let mut writer = bam::io::Writer::new(io::sink());
    /// writer.write_header(&header)?;
    ///
    /// let record = bam::Record::default();
    /// writer.write_record(&header, &record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, header: &sam::Header, record: &Record) -> io::Result<()> {
        self.write_alignment_record(header, record)
    }
}

impl<W> Writer<bgzf::io::Writer<W>>
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
    /// # use std::io;
    /// use noodles_bam as bam;
    /// let writer = bam::io::Writer::new(io::sink());
    /// ```
    pub fn new(writer: W) -> Self {
        Self::from(bgzf::io::Writer::new(writer))
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
    /// let mut writer = bam::io::Writer::new(io::sink());
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
        use crate::record::codec::encode;

        self.buf.clear();
        encode(&mut self.buf, header, record)?;

        let block_size = u32::try_from(self.buf.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_u32_le(&mut self.inner, block_size)?;

        self.inner.write_all(&self.buf)?;

        Ok(())
    }

    fn finish(&mut self, _: &sam::Header) -> io::Result<()> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use sam::alignment::{
        RecordBuf,
        record::Flags,
        record_buf::{QualityScores, Sequence},
    };

    use super::*;
    use crate::io::Reader;

    #[test]
    fn test_write_alignment_record() -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let record = RecordBuf::default();
        writer.write_alignment_record(&header, &record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().get_ref().as_slice());

        let mut record = RecordBuf::default();
        reader.read_record_buf(&header, &mut record)?;

        assert!(record.name().is_none());
        assert_eq!(record.flags(), Flags::UNMAPPED);
        assert!(record.reference_sequence_id().is_none());
        assert!(record.alignment_start().is_none());
        assert!(record.mapping_quality().is_none());
        assert!(record.cigar().as_ref().is_empty());
        assert!(record.mate_reference_sequence_id().is_none());
        assert!(record.mate_alignment_start().is_none());
        assert_eq!(record.template_length(), 0);
        assert!(record.sequence().is_empty());
        assert!(record.quality_scores().is_empty());
        assert!(record.data().is_empty());

        Ok(())
    }

    #[test]
    fn test_write_alignment_record_with_sequence_length_less_than_quality_scores_length()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();

        let mut record = RecordBuf::builder()
            .set_sequence(Sequence::from(b"AT"))
            .build();

        *record.quality_scores_mut() = QualityScores::from(vec![45, 35, 43, 50]);

        assert!(writer.write_alignment_record(&header, &record).is_err());

        Ok(())
    }

    #[test]
    fn test_write_alignment_record_with_sequence_length_greater_than_quality_scores_length()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();

        let mut record = RecordBuf::builder()
            .set_sequence(Sequence::from(b"ATCG"))
            .build();

        *record.quality_scores_mut() = QualityScores::from(vec![45, 35]);

        assert!(writer.write_alignment_record(&header, &record).is_err());

        Ok(())
    }

    #[test]
    fn test_write_alignment_record_with_no_sequence_and_with_quality_scores()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let mut record = RecordBuf::default();
        *record.quality_scores_mut() = QualityScores::from(vec![45, 35, 43, 50]);

        assert!(writer.write_alignment_record(&header, &record).is_err());

        Ok(())
    }

    #[test]
    fn test_write_alignment_record_with_sequence_and_no_quality_scores()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();

        let record = RecordBuf::builder()
            .set_sequence(Sequence::from(b"ATCG"))
            .build();

        writer.write_alignment_record(&header, &record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().get_ref().as_slice());

        let mut record = RecordBuf::default();
        reader.read_record_buf(&header, &mut record)?;

        let expected = Sequence::from(b"ATCG");
        assert_eq!(record.sequence(), &expected);

        assert!(record.quality_scores().is_empty());

        Ok(())
    }

    #[test]
    fn test_write_alignment_record_with_sequence_and_quality_scores()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let sam_record = RecordBuf::builder()
            .set_sequence(Sequence::from(b"ATCG"))
            .set_quality_scores(QualityScores::from(vec![45, 35, 43, 50]))
            .build();

        writer.write_alignment_record(&header, &sam_record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().get_ref().as_slice());

        let mut record = RecordBuf::default();
        reader.read_record_buf(&header, &mut record)?;

        let expected = Sequence::from(b"ATCG");
        assert_eq!(record.sequence(), &expected);

        assert_eq!(record.quality_scores(), sam_record.quality_scores());

        Ok(())
    }

    #[test]
    fn test_write_alignment_record_with_data() -> io::Result<()> {
        use sam::alignment::{record::data::field::Tag, record_buf::data::field::Value};

        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let sam_record = RecordBuf::builder()
            .set_data(
                [
                    (Tag::READ_GROUP, Value::from("rg0")),
                    (Tag::ALIGNMENT_HIT_COUNT, Value::UInt8(1)),
                ]
                .into_iter()
                .collect(),
            )
            .build();

        writer.write_alignment_record(&header, &sam_record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().get_ref().as_slice());

        let mut record = RecordBuf::default();
        reader.read_record_buf(&header, &mut record)?;

        let bam_data = record.data();
        let mut fields = bam_data.iter();

        assert_eq!(fields.next(), Some((Tag::READ_GROUP, &Value::from("rg0"))));

        assert_eq!(
            fields.next(),
            Some((Tag::ALIGNMENT_HIT_COUNT, &Value::UInt8(1)))
        );

        assert!(fields.next().is_none());

        Ok(())
    }
}
