use std::{mem, ops::{Deref, DerefMut}};

use byteorder::{ByteOrder, LittleEndian};

use super::{Cigar, Data, Flag, Quality, Sequence};

#[derive(Clone, Default)]
pub struct Record(Vec<u8>);

impl Record {
    pub fn new() -> Record {
        Record::default()
    }

    pub fn with_capacity(capacity: usize) -> Record {
        Record(Vec::with_capacity(capacity))
    }

    pub fn resize(&mut self, new_len: usize) {
        self.0.resize(new_len, Default::default());
    }

    pub fn block_size(&self) -> i32 {
        self.0.len() as i32
    }

    pub fn ref_id(&self) -> i32 {
        let offset = 0;
        let len = mem::size_of::<i32>();
        LittleEndian::read_i32(&self.0[offset..offset + len])
    }

    pub fn pos(&self) -> i32 {
        let offset = 4;
        let len = mem::size_of::<i32>();
        LittleEndian::read_i32(&self.0[offset..offset + len])
    }

    pub fn l_read_name(&self) -> u8 {
        let offset = 8;
        self.0[offset]
    }

    pub fn mapq(&self) -> u8 {
        let offset = 9;
        self.0[offset]
    }

    pub fn bin(&self) -> u16 {
        let offset = 10;
        let len = mem::size_of::<u16>();
        LittleEndian::read_u16(&self.0[offset..offset + len])
    }

    pub fn n_cigar_op(&self) -> u16 {
        let offset = 12;
        let len = mem::size_of::<u16>();
        LittleEndian::read_u16(&self.0[offset..offset + len])
    }

    pub fn flag(&self) -> Flag {
        let offset = 14;
        let len = mem::size_of::<u16>();
        let value = LittleEndian::read_u16(&self.0[offset..offset + len]);
        Flag::from(value)
    }

    pub fn l_seq(&self) -> i32 {
        let offset = 16;
        let len = mem::size_of::<i32>();
        LittleEndian::read_i32(&self.0[offset..offset + len])
    }

    pub fn next_ref_id(&self) -> i32 {
        let offset = 20;
        let len = mem::size_of::<i32>();
        LittleEndian::read_i32(&self.0[offset..offset + len])
    }

    pub fn next_pos(&self) -> i32 {
        let offset = 24;
        let len = mem::size_of::<i32>();
        LittleEndian::read_i32(&self.0[offset..offset + len])
    }

    pub fn tlen(&self) -> i32 {
        let offset = 28;
        let len = mem::size_of::<i32>();
        LittleEndian::read_i32(&self.0[offset..offset + len])
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

    pub fn seq(&self) -> Sequence {
        let offset = 32
            + (self.l_read_name() as usize)
            + mem::size_of::<u32>() * (self.n_cigar_op() as usize);
        let len = ((self.l_seq() + 1) / 2) as usize;

        let bytes = &self.0[offset..offset + len];
        let n_chars = self.l_seq() as usize;
        Sequence::new(bytes, n_chars)
    }

    pub fn qual(&self) -> Quality {
        let l_seq = self.l_seq();

        let offset = 32
            + (self.l_read_name() as usize)
            + mem::size_of::<u32>() * (self.n_cigar_op() as usize)
            + ((self.l_seq() + 1) / 2) as usize;
        let len = l_seq as usize;

        let bytes = &self.0[offset..offset + len];
        Quality::new(bytes)
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

impl From<Vec<u8>> for Record {
    fn from(bytes: Vec<u8>) -> Record {
        Record(bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::Record;

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
        let r = Record::new();
        assert_eq!(r.block_size(), 0);

        let r = build_record();
        assert_eq!(r.block_size(), 64);
    }

    #[test]
    fn test_ref_id() {
        let r = build_record();
        assert_eq!(r.ref_id(), 10);
    }

    #[test]
    fn test_pos() {
        let r = build_record();
        assert_eq!(r.pos(), 61061);
    }

    #[test]
    fn test_l_read_name() {
        let r = build_record();
        assert_eq!(r.l_read_name(), 10);
    }

    #[test]
    fn test_mapq() {
        let r = build_record();
        assert_eq!(r.mapq(), 12);
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
    fn test_flag() {
        let r = build_record();
        assert_eq!(r.flag().inner(), 0xa3);
    }

    #[test]
    fn test_l_seq() {
        let r = build_record();
        assert_eq!(r.l_seq(), 4);
    }

    #[test]
    fn test_next_ref_id() {
        let r = build_record();
        assert_eq!(r.next_ref_id(), 10);
    }

    #[test]
    fn test_next_pos() {
        let r = build_record();
        assert_eq!(r.next_pos(), 61152);
    }

    #[test]
    fn test_tlen() {
        let r = build_record();
        assert_eq!(r.tlen(), 166);
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
    fn test_seq() {
        let r = build_record();
        let expected = [0x18, 0x42];
        assert_eq!(*r.seq(), expected);
    }

    #[test]
    fn test_qual() {
        let r = build_record();
        let expected = [0x1f, 0x1d, 0x1e, 0x20];
        assert_eq!(*r.qual(), expected);
    }

    #[test]
    fn test_data() {
        let r = build_record();
        let expected = [0x4e, 0x4d, 0x43, 0x00, 0x50, 0x47, 0x5a, 0x53, 0x4e, 0x41, 0x50, 0x00];
        assert_eq!(*r.data(), expected);
    }
}
