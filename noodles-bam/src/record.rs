pub mod cigar;
pub mod data;
mod quality_scores;
pub mod sequence;

pub use self::{cigar::Cigar, data::Data, quality_scores::QualityScores, sequence::Sequence};

use std::{
    ffi::CStr,
    fmt, mem,
    ops::{Deref, DerefMut},
};

use byteorder::{ByteOrder, LittleEndian};
use noodles_sam as sam;

#[derive(Clone)]
pub struct Record(Vec<u8>);

impl Record {
    pub fn with_capacity(capacity: usize) -> Self {
        Self(Vec::with_capacity(capacity))
    }

    pub fn resize(&mut self, new_len: usize) {
        self.0.resize(new_len, Default::default());
    }

    pub fn block_size(&self) -> u32 {
        self.0.len() as u32
    }

    pub fn reference_sequence_id(&self) -> i32 {
        LittleEndian::read_i32(&self.0)
    }

    pub fn position(&self) -> i32 {
        let offset = 4;
        LittleEndian::read_i32(&self.0[offset..])
    }

    fn l_read_name(&self) -> u8 {
        let offset = 8;
        self.0[offset]
    }

    pub fn mapping_quality(&self) -> sam::record::MappingQuality {
        let offset = 9;
        sam::record::MappingQuality::from(self.0[offset])
    }

    pub fn bin(&self) -> u16 {
        let offset = 10;
        LittleEndian::read_u16(&self.0[offset..])
    }

    fn n_cigar_op(&self) -> u16 {
        let offset = 12;
        LittleEndian::read_u16(&self.0[offset..])
    }

    pub fn flags(&self) -> sam::record::Flags {
        let offset = 14;
        let value = LittleEndian::read_u16(&self.0[offset..]);
        sam::record::Flags::from(value)
    }

    fn l_seq(&self) -> u32 {
        let offset = 16;
        LittleEndian::read_u32(&self.0[offset..])
    }

    pub fn mate_reference_sequence_id(&self) -> i32 {
        let offset = 20;
        LittleEndian::read_i32(&self.0[offset..])
    }

    pub fn mate_position(&self) -> i32 {
        let offset = 24;
        LittleEndian::read_i32(&self.0[offset..])
    }

    pub fn template_len(&self) -> i32 {
        let offset = 28;
        LittleEndian::read_i32(&self.0[offset..])
    }

    pub fn read_name(&self) -> &[u8] {
        let offset = 32;
        let len = self.l_read_name() as usize;
        &self.0[offset..offset + len]
    }

    pub fn cigar(&self) -> Cigar {
        let offset = 32 + (self.l_read_name() as usize);
        let len = mem::size_of::<u32>() * (self.n_cigar_op() as usize);
        let bytes = &self.0[offset..offset + len];
        Cigar::new(bytes)
    }

    pub fn sequence(&self) -> Sequence {
        let offset = 32
            + (self.l_read_name() as usize)
            + mem::size_of::<u32>() * (self.n_cigar_op() as usize);
        let len = ((self.l_seq() + 1) / 2) as usize;

        let bytes = &self.0[offset..offset + len];
        let n_chars = self.l_seq() as usize;
        Sequence::new(bytes, n_chars)
    }

    pub fn quality_scores(&self) -> QualityScores {
        let l_seq = self.l_seq();

        let offset = 32
            + (self.l_read_name() as usize)
            + mem::size_of::<u32>() * (self.n_cigar_op() as usize)
            + ((self.l_seq() + 1) / 2) as usize;
        let len = l_seq as usize;

        let bytes = &self.0[offset..offset + len];
        QualityScores::new(bytes)
    }

    pub fn data(&self) -> Data {
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
            0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x02, 0xff, 0x48, 0x12, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
            0x00, 0x00, 0x00, 0x00, 0x2a, 0x00,
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
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        let read_name = CStr::from_bytes_with_nul(self.read_name());

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
            .field("tlen", &self.template_len())
            .field("read_name", &read_name)
            .field("cigar", &self.cigar())
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
    use super::*;

    #[rustfmt::skip]
    fn build_record() -> Record {
        Record::from(vec![
            // ref_id
            0x0a, 0x00, 0x00, 0x00,
            // pos
            0x85, 0xee, 0x00, 0x00,
            // l_read_name
            0x0a,
            // mapq
            0x0c,
            // bin
            0x4c, 0x12,
            // n_ciar_op
            0x01, 0x00,
            // flag
            0xa3, 0x00,
            // l_seq
            0x04, 0x00, 0x00, 0x00,
            // next_ref_id
            0x0a, 0x00, 0x00, 0x00,
            // next_pos
            0xe0, 0xee, 0x00, 0x00,
            // tlen
            0xa6, 0x00, 0x00, 0x00,
            // read_name
            0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73, 0x3a, 0x30, 0x00,
            // cigar
            0x40, 0x00, 0x00, 0x00,
            // seq
            0x18, 0x42,
            // qual
            0x1f, 0x1d, 0x1e, 0x20,
            // data
            0x4e, 0x4d, 0x43, 0x00, 0x50, 0x47, 0x5a, 0x53, 0x4e, 0x41, 0x50, 0x00,
        ])
    }

    #[test]
    fn test_block_size() {
        let r = build_record();
        assert_eq!(r.block_size(), 64);
    }

    #[test]
    fn test_reference_sequence_id() {
        let r = build_record();
        assert_eq!(r.reference_sequence_id(), 10);
    }

    #[test]
    fn test_position() {
        let r = build_record();
        assert_eq!(r.position(), 61061);
    }

    #[test]
    fn test_l_read_name() {
        let r = build_record();
        assert_eq!(r.l_read_name(), 10);
    }

    #[test]
    fn test_mapping_quality() {
        let r = build_record();
        assert_eq!(r.mapping_quality(), sam::record::MappingQuality::from(12));
    }

    #[test]
    fn test_bin() {
        let r = build_record();
        assert_eq!(r.bin(), 4684);
    }

    #[test]
    fn test_n_cigar_op() {
        let r = build_record();
        assert_eq!(r.n_cigar_op(), 1);
    }

    #[test]
    fn test_flags() {
        let r = build_record();
        assert_eq!(u16::from(r.flags()), 0xa3);
    }

    #[test]
    fn test_l_seq() {
        let r = build_record();
        assert_eq!(r.l_seq(), 4);
    }

    #[test]
    fn test_mate_reference_sequence_id() {
        let r = build_record();
        assert_eq!(r.mate_reference_sequence_id(), 10);
    }

    #[test]
    fn test_mate_position() {
        let r = build_record();
        assert_eq!(r.mate_position(), 61152);
    }

    #[test]
    fn test_template_len() {
        let r = build_record();
        assert_eq!(r.template_len(), 166);
    }

    #[test]
    fn test_read_name() {
        let r = build_record();
        let expected = [0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73, 0x3a, 0x30, 0x00];
        assert_eq!(r.read_name(), expected);
    }

    #[test]
    fn test_cigar() {
        let r = build_record();
        let expected = [0x40, 0x00, 0x00, 0x00];
        assert_eq!(*r.cigar(), expected);
    }

    #[test]
    fn test_sequence() {
        let r = build_record();
        let expected = [0x18, 0x42];
        assert_eq!(*r.sequence(), expected);
    }

    #[test]
    fn test_quality_scores() {
        let r = build_record();
        let expected = [0x1f, 0x1d, 0x1e, 0x20];
        assert_eq!(*r.quality_scores(), expected);
    }

    #[test]
    fn test_data() {
        let r = build_record();
        let expected = [
            0x4e, 0x4d, 0x43, 0x00, 0x50, 0x47, 0x5a, 0x53, 0x4e, 0x41, 0x50, 0x00,
        ];
        assert_eq!(*r.data(), expected);
    }
}
