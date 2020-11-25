//! BAM record and fields.

pub mod cigar;
mod convert;
pub mod data;
pub mod quality_scores;
pub mod reference_sequence_id;
pub mod sequence;

pub use self::{
    cigar::Cigar, data::Data, quality_scores::QualityScores,
    reference_sequence_id::ReferenceSequenceId, sequence::Sequence,
};

use std::{
    convert::TryFrom,
    ffi::{self, CStr},
    fmt, mem,
    ops::{Deref, DerefMut},
};

use byteorder::{ByteOrder, LittleEndian};
use noodles_sam as sam;

pub(crate) const UNMAPPED_POSITION: i32 = -1;

/// A BAM record.
///
/// A BAM record encodes the same fields as a SAM record:
///
///  * reference sequence ID (`RNAME` equiv.),
///  * position (`POS`),
///  * mapping quality (`MAPQ`),
///  * flags (`FLAG`),
///  * mate reference sequence ID (`RNEXT` equiv.),
///  * mate position (`PNEXT`),
///  * template length (`TLEN`),
///  * read name (`QNAME`),
///  * CIGAR operations (`CIGAR`),
///  * sequence (`SEQ`),
///  * quality scores (`QUAL`), and
///  * optional data fields.
///
/// Additionally, it encodes the BAM index bin (`bin`).
///
/// A `bam::Record` wraps a raw byte buffer, and the fields should be considered immutable.
#[derive(Clone, Eq, PartialEq)]
pub struct Record(Vec<u8>);

impl Record {
    pub(crate) fn resize(&mut self, new_len: usize) {
        self.0.resize(new_len, Default::default());
    }

    pub(crate) fn block_size(&self) -> u32 {
        self.0.len() as u32
    }

    /// Returns the reference sequence ID of this record.
    ///
    /// The reference sequence ID is the index of the associated reference sequence in the SAM
    /// header or BAM reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.reference_sequence_id().is_none());
    /// ```
    pub fn reference_sequence_id(&self) -> Option<ReferenceSequenceId> {
        let id = LittleEndian::read_i32(&self.0);

        if id == reference_sequence_id::UNMAPPED {
            None
        } else {
            ReferenceSequenceId::try_from(id).ok()
        }
    }

    /// Returns the start position of this record.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.position().is_none());
    /// ```
    pub fn position(&self) -> Option<sam::record::Position> {
        let offset = 4;
        let pos = LittleEndian::read_i32(&self.0[offset..]);

        if pos == UNMAPPED_POSITION {
            None
        } else {
            sam::record::Position::try_from(pos + 1).ok()
        }
    }

    fn l_read_name(&self) -> u8 {
        let offset = 8;
        self.0[offset]
    }

    /// Returns the mapping quality of this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    /// let record = bam::Record::default();
    /// assert!(record.mapping_quality().is_none());
    /// ```
    pub fn mapping_quality(&self) -> sam::record::MappingQuality {
        let offset = 9;
        sam::record::MappingQuality::from(self.0[offset])
    }

    /// Returns the index bin that includes this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert_eq!(record.bin(), 4680);
    /// ```
    pub fn bin(&self) -> u16 {
        let offset = 10;
        LittleEndian::read_u16(&self.0[offset..])
    }

    fn n_cigar_op(&self) -> u16 {
        let offset = 12;
        LittleEndian::read_u16(&self.0[offset..])
    }

    /// Returns the SAM flags of this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    /// let record = bam::Record::default();
    /// assert_eq!(record.flags(), sam::record::Flags::UNMAPPED);
    /// ```
    pub fn flags(&self) -> sam::record::Flags {
        let offset = 14;
        let value = LittleEndian::read_u16(&self.0[offset..]);
        sam::record::Flags::from(value)
    }

    fn l_seq(&self) -> u32 {
        let offset = 16;
        LittleEndian::read_u32(&self.0[offset..])
    }

    /// Returns the reference sequence ID of the mate of this record.
    ///
    /// The mate reference sequence ID is the index of the associated reference sequence in the SAM
    /// header or BAM reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.mate_reference_sequence_id().is_none());
    /// ```
    pub fn mate_reference_sequence_id(&self) -> Option<ReferenceSequenceId> {
        let offset = 20;
        let id = LittleEndian::read_i32(&self.0[offset..]);

        if id == reference_sequence_id::UNMAPPED {
            None
        } else {
            ReferenceSequenceId::try_from(id).ok()
        }
    }

    /// Returns the start position of the mate of this record.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.mate_position().is_none());
    /// ```
    pub fn mate_position(&self) -> Option<sam::record::Position> {
        let offset = 24;
        let pos = LittleEndian::read_i32(&self.0[offset..]);

        if pos == UNMAPPED_POSITION {
            None
        } else {
            sam::record::Position::try_from(pos + 1).ok()
        }
    }

    /// Returns the template length of this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert_eq!(record.template_length(), 0);
    /// ```
    pub fn template_length(&self) -> i32 {
        let offset = 28;
        LittleEndian::read_i32(&self.0[offset..])
    }

    /// Returns the read name of this record.
    ///
    /// This is also called the query name.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::ffi;
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert_eq!(record.read_name()?.to_bytes(), b"*");
    /// # Ok::<(), ffi::FromBytesWithNulError>(())
    /// ```
    pub fn read_name(&self) -> Result<&CStr, ffi::FromBytesWithNulError> {
        let offset = 32;
        let len = self.l_read_name() as usize;
        let data = &self.0[offset..offset + len];
        CStr::from_bytes_with_nul(data)
    }

    /// Returns the CIGAR operations that describe how the read was mapped.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.cigar().is_empty());
    /// ```
    pub fn cigar(&self) -> Cigar<'_> {
        let offset = 32 + (self.l_read_name() as usize);
        let len = mem::size_of::<u32>() * (self.n_cigar_op() as usize);
        let bytes = &self.0[offset..offset + len];
        Cigar::new(bytes)
    }

    /// Returns the bases in the sequence of this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.sequence().is_empty());
    /// ```
    pub fn sequence(&self) -> Sequence<'_> {
        let offset = 32
            + (self.l_read_name() as usize)
            + mem::size_of::<u32>() * (self.n_cigar_op() as usize);
        let len = ((self.l_seq() + 1) / 2) as usize;

        let bytes = &self.0[offset..offset + len];
        let base_count = self.l_seq() as usize;
        Sequence::new(bytes, base_count)
    }

    /// Returns the quality score for each base in the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.quality_scores().is_empty());
    /// ```
    pub fn quality_scores(&self) -> QualityScores<'_> {
        let l_seq = self.l_seq();

        let offset = 32
            + (self.l_read_name() as usize)
            + mem::size_of::<u32>() * (self.n_cigar_op() as usize)
            + ((self.l_seq() + 1) / 2) as usize;
        let len = l_seq as usize;

        let bytes = &self.0[offset..offset + len];
        QualityScores::new(bytes)
    }

    /// Returns the optional data fields for this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.data().is_empty());
    /// ```
    pub fn data(&self) -> Data<'_> {
        let l_seq = self.l_seq();

        let offset = 32
            + (self.l_read_name() as usize)
            + mem::size_of::<u32>() * (self.n_cigar_op() as usize)
            + ((self.l_seq() + 1) / 2) as usize
            + l_seq as usize;
        let len = self.block_size() as usize;

        let bytes = &self.0[offset..len];
        Data::new(bytes)
    }
}

impl Default for Record {
    fn default() -> Self {
        Self::from(vec![
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
        ])
    }
}

impl Deref for Record {
    type Target = [u8];

    fn deref(&self) -> &[u8] {
        &self.0
    }
}

impl DerefMut for Record {
    fn deref_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl fmt::Debug for Record {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.debug_struct("Record")
            .field("block_size", &self.block_size())
            .field("ref_id", &self.reference_sequence_id())
            .field("pos", &self.position())
            .field("l_read_name", &self.l_read_name())
            .field("mapq", &self.mapping_quality())
            .field("bin", &self.bin())
            .field("n_cigar_op", &self.n_cigar_op())
            .field("flag", &self.flags())
            .field("l_seq", &self.l_seq())
            .field("next_ref_id", &self.mate_reference_sequence_id())
            .field("next_pos", &self.mate_position())
            .field("tlen", &self.template_length())
            .field("read_name", &self.read_name())
            .field("cigar", &self.cigar())
            .field("seq", &self.sequence())
            .finish()
    }
}

impl From<Vec<u8>> for Record {
    fn from(bytes: Vec<u8>) -> Self {
        Self(bytes)
    }
}

#[cfg(test)]
mod tests {
    use std::{
        ffi::CString,
        io::{self, BufWriter, Write},
    };

    use byteorder::WriteBytesExt;

    use super::*;

    fn build_record() -> io::Result<Record> {
        use sam::record::Flags;

        let read_name = CString::new("noodles:0")
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let flag = u16::from(Flags::PAIRED | Flags::READ_1);

        let mut writer = BufWriter::new(Vec::new());

        // ref_id
        writer.write_i32::<LittleEndian>(10)?;
        // pos
        writer.write_i32::<LittleEndian>(61061)?;
        // l_read_name
        writer.write_u8(read_name.as_bytes_with_nul().len() as u8)?;
        // mapq
        writer.write_u8(12)?;
        // bin
        writer.write_u16::<LittleEndian>(4684)?;
        // n_ciar_op
        writer.write_u16::<LittleEndian>(1)?;
        // flag
        writer.write_u16::<LittleEndian>(flag)?;
        // l_seq
        writer.write_u32::<LittleEndian>(4)?;
        // next_ref_id
        writer.write_i32::<LittleEndian>(10)?;
        // next_pos
        writer.write_i32::<LittleEndian>(61152)?;
        // tlen
        writer.write_i32::<LittleEndian>(166)?;
        // read_name
        writer.write_all(read_name.as_bytes_with_nul())?;
        // cigar
        writer.write_all(&[0x40, 0x00, 0x00, 0x00])?;
        // seq
        writer.write_all(&[0x18, 0x42])?;
        // qual
        writer.write_all(&[0x1f, 0x1d, 0x1e, 0x20])?;
        // data
        writer.write_all(&[
            0x4e, 0x4d, 0x43, 0x00, 0x50, 0x47, 0x5a, 0x53, 0x4e, 0x41, 0x50, 0x00,
        ])?;

        writer
            .into_inner()
            .map(Record::from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    }

    #[test]
    fn test_block_size() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.block_size(), 64);
        Ok(())
    }

    #[test]
    fn test_reference_sequence_id() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.reference_sequence_id().map(i32::from), Some(10));
        Ok(())
    }

    #[test]
    fn test_position() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.position().map(i32::from), Some(61062));
        Ok(())
    }

    #[test]
    fn test_l_read_name() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.l_read_name(), 10);
        Ok(())
    }

    #[test]
    fn test_mapping_quality() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(
            record.mapping_quality(),
            sam::record::MappingQuality::from(12)
        );
        Ok(())
    }

    #[test]
    fn test_bin() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.bin(), 4684);
        Ok(())
    }

    #[test]
    fn test_n_cigar_op() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.n_cigar_op(), 1);
        Ok(())
    }

    #[test]
    fn test_flags() -> io::Result<()> {
        use sam::record::Flags;
        let record = build_record()?;
        assert_eq!(record.flags(), Flags::PAIRED | Flags::READ_1);
        Ok(())
    }

    #[test]
    fn test_l_seq() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.l_seq(), 4);
        Ok(())
    }

    #[test]
    fn test_mate_reference_sequence_id() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.mate_reference_sequence_id().map(i32::from), Some(10));
        Ok(())
    }

    #[test]
    fn test_mate_position() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.mate_position().map(i32::from), Some(61153));
        Ok(())
    }

    #[test]
    fn test_template_length() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.template_length(), 166);
        Ok(())
    }

    #[test]
    fn test_read_name() -> Result<(), Box<dyn std::error::Error>> {
        let record = build_record()?;
        assert_eq!(record.read_name()?.to_bytes(), b"noodles:0");
        Ok(())
    }

    #[test]
    fn test_cigar() -> io::Result<()> {
        let record = build_record()?;
        let expected = [0x40, 0x00, 0x00, 0x00];
        assert_eq!(*record.cigar(), expected);
        Ok(())
    }

    #[test]
    fn test_sequence() -> io::Result<()> {
        let record = build_record()?;
        let expected = [0x18, 0x42];
        assert_eq!(*record.sequence(), expected);
        Ok(())
    }

    #[test]
    fn test_quality_scores() -> io::Result<()> {
        let record = build_record()?;
        let expected = [0x1f, 0x1d, 0x1e, 0x20];
        assert_eq!(*record.quality_scores(), expected);
        Ok(())
    }

    #[test]
    fn test_data() -> io::Result<()> {
        let record = build_record()?;
        let expected = [
            0x4e, 0x4d, 0x43, 0x00, 0x50, 0x47, 0x5a, 0x53, 0x4e, 0x41, 0x50, 0x00,
        ];
        assert_eq!(*record.data(), expected);
        Ok(())
    }
}
