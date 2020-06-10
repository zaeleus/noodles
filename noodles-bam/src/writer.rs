use std::{
    ffi::CString,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_bgzf as bgzf;
use noodles_sam::{
    self as sam,
    header::{ReferenceSequence, ReferenceSequences},
    record::{Cigar, MateReferenceSequenceName, QualityScores, Sequence},
};

use super::MAGIC_NUMBER;

// § 4.2 The BAM format (2020-04-30)
//
// ref_id (4) + pos (4) + l_read_name (1) + mapq (1) + bin (2) + n_cigar_op (2) + flag (2) + l_seq
// (4) + next_ref_id (4) + next_pos (4) + tlen (4)
const BLOCK_HEADER_SIZE: usize = 32;

// § 4.2.1 BIN field calculation (2020-04-30)
const UNMAPPED_BIN: u16 = 4680;

// § 4.2.3 SEQ and QUAL encoding (2020-04-30)
const NULL_QUALITY_SCORE: u8 = 255;

/// A BAM writer.
///
/// Since the raw text header and `bam::Record` are immutable, BAM files are created by encoding a
/// SAM header and SAM records.
///
/// # Examples
///
/// ```no_run
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
/// writer.write_record(header.reference_sequences(), &record)?;
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
    /// writer.write_record(&reference_sequences, &record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(
        &mut self,
        reference_sequences: &ReferenceSequences,
        record: &sam::Record,
    ) -> io::Result<()> {
        let c_read_name = CString::new(record.read_name().as_ref())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let reference_sequence_id = match &**record.reference_sequence_name() {
            Some(name) => reference_sequences
                .get_full(name)
                .map(|(i, _, _)| i as i32)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "invalid reference sequence id")
                })?,
            None => -1,
        };

        let mate_reference_sequence_id = match record.mate_reference_sequence_name() {
            MateReferenceSequenceName::Some(name) => reference_sequences
                .get_full(name)
                .map(|(i, _, _)| i as i32)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "invalid reference sequence id")
                })?,
            MateReferenceSequenceName::Eq => reference_sequence_id,
            MateReferenceSequenceName::None => -1,
        };

        let read_name = c_read_name.as_bytes_with_nul();
        let l_read_name = read_name.len() as u8;
        let n_cigar_op = record.cigar().ops().len() as u16;
        let l_seq = record.sequence().len() as i32;

        let block_size = BLOCK_HEADER_SIZE as i32
            + (l_read_name as i32)
            + (4 * (n_cigar_op as i32))
            + ((l_seq + 1) / 2)
            + l_seq;

        self.inner.write_i32::<LittleEndian>(block_size)?;

        let ref_id = reference_sequence_id as i32;
        self.inner.write_i32::<LittleEndian>(ref_id)?;

        let pos = i32::from(record.position()) - 1;
        self.inner.write_i32::<LittleEndian>(pos)?;

        self.inner.write_u8(l_read_name)?;

        let mapq = u8::from(record.mapping_quality());
        self.inner.write_u8(mapq)?;

        let bin = record
            .position()
            .map(|start| {
                let end = record.cigar().mapped_len() as i32;
                region_to_bin(start, end) as u16
            })
            .unwrap_or(UNMAPPED_BIN);

        self.inner.write_u16::<LittleEndian>(bin)?;

        self.inner.write_u16::<LittleEndian>(n_cigar_op)?;

        let flag = u16::from(record.flags());
        self.inner.write_u16::<LittleEndian>(flag)?;

        self.inner.write_i32::<LittleEndian>(l_seq)?;

        let next_ref_id = mate_reference_sequence_id as i32;
        self.inner.write_i32::<LittleEndian>(next_ref_id)?;

        let next_pos = i32::from(record.mate_position()) - 1;
        self.inner.write_i32::<LittleEndian>(next_pos)?;

        let tlen = record.template_len();
        self.inner.write_i32::<LittleEndian>(tlen)?;

        self.inner.write_all(read_name)?;

        write_cigar(&mut self.inner, record.cigar())?;

        // § 4.2.3 SEQ and QUAL encoding (2020-04-30)
        let sequence = record.sequence();
        let quality_scores = record.quality_scores();

        write_seq(&mut self.inner, sequence)?;

        if sequence.len() < quality_scores.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "quality scores length does not match sequence length",
            ));
        } else if sequence.len() > quality_scores.len() {
            if quality_scores.is_empty() {
                for _ in 0..sequence.len() {
                    self.inner.write_u8(NULL_QUALITY_SCORE)?;
                }
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "quality scores length does not match sequence length",
                ));
            }
        } else {
            write_qual(&mut self.inner, quality_scores)?;
        }

        Ok(())
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

fn write_cigar<W>(writer: &mut W, cigar: &Cigar) -> io::Result<()>
where
    W: Write,
{
    for op in cigar.ops() {
        let len = op.len() as u32;
        let kind = op.kind() as u32;
        let value = len << 4 | kind;
        writer.write_u32::<LittleEndian>(value)?;
    }

    Ok(())
}

fn write_seq<W>(writer: &mut W, sequence: &Sequence) -> io::Result<()>
where
    W: Write,
{
    use sam::record::sequence::Base as SamBase;

    fn base_to_u8(base: SamBase) -> u8 {
        match base {
            SamBase::Eq => 0,
            SamBase::A => 1,
            SamBase::C => 2,
            SamBase::M => 3,
            SamBase::G => 4,
            SamBase::R => 5,
            SamBase::S => 6,
            SamBase::V => 7,
            SamBase::T => 8,
            SamBase::W => 9,
            SamBase::Y => 10,
            SamBase::H => 11,
            SamBase::K => 12,
            SamBase::D => 13,
            SamBase::B => 14,
            _ => 15,
        }
    }

    for chunk in sequence.bases().chunks(2) {
        let l = base_to_u8(chunk[0]);

        let r = if let Some(c) = chunk.get(1) {
            base_to_u8(*c)
        } else {
            0
        };

        let value = l << 4 | r;

        writer.write_u8(value)?;
    }

    Ok(())
}

fn write_qual<W>(writer: &mut W, quality_scores: &QualityScores) -> io::Result<()>
where
    W: Write,
{
    for score in quality_scores.scores() {
        let value = u8::from(*score);
        writer.write_u8(value)?;
    }

    Ok(())
}

// See § 5.3 in SAMv1.pdf (accessed 2020-04-24).
#[allow(clippy::eq_op)]
fn region_to_bin(start: i32, mut end: i32) -> i32 {
    end -= 1;

    if start >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (start >> 14)
    } else if start >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (start >> 17)
    } else if start >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (start >> 20)
    } else if start >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (start >> 23)
    } else if start >> 26 == end >> 26 {
        ((1 << 3) - 1) / 7 + (start >> 26)
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
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
    fn test_write_record() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let sam_record = sam::Record::default();
        writer.write_record(header.reference_sequences(), &sam_record)?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().as_slice());

        let mut record = Record::default();
        reader.read_record(&mut record)?;

        assert_eq!(record.read_name(), b"*\0");
        assert_eq!(record.flags(), sam::record::Flags::default());
        assert_eq!(record.reference_sequence_id(), -1);
        assert_eq!(record.position(), -1);
        assert_eq!(
            record.mapping_quality(),
            sam::record::MappingQuality::from(255)
        );
        assert!(record.cigar().is_empty());
        assert_eq!(record.mate_reference_sequence_id(), -1);
        assert_eq!(record.mate_position(), -1);
        assert_eq!(record.template_len(), 0);
        assert!(record.sequence().is_empty());
        assert!(record.quality_scores().is_empty());
        assert!(record.data().is_empty());

        Ok(())
    }

    #[test]
    fn test_write_record_with_sequence_length_less_than_quality_scores_length(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let record = sam::Record::builder()
            .set_sequence("AT".parse()?)
            .set_quality_scores("NDLS".parse()?)
            .build();

        assert!(writer
            .write_record(header.reference_sequences(), &record)
            .is_err());

        Ok(())
    }

    #[test]
    fn test_write_record_with_sequence_length_greater_than_quality_scores_length(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let record = sam::Record::builder()
            .set_sequence("ATCG".parse()?)
            .set_quality_scores("ND".parse()?)
            .build();

        assert!(writer
            .write_record(header.reference_sequences(), &record)
            .is_err());

        Ok(())
    }

    #[test]
    fn test_write_record_with_sequence_and_no_quality_scores(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let sam_record = sam::Record::builder().set_sequence("ATCG".parse()?).build();
        writer.write_record(header.reference_sequences(), &sam_record)?;

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
    fn test_write_record_with_sequence_and_quality_scores() -> Result<(), Box<dyn std::error::Error>>
    {
        let mut writer = Writer::new(Vec::new());

        let header = sam::Header::default();
        let sam_record = sam::Record::builder()
            .set_sequence("ATCG".parse()?)
            .set_quality_scores("NDLS".parse()?)
            .build();

        writer.write_record(header.reference_sequences(), &sam_record)?;
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
}
