use formats::bam::{Cigar, Data, Quality, Sequence};

#[derive(Debug)]
pub struct Flag(u16);

impl Flag {
    pub fn new(flag: u16) -> Flag {
        Flag(flag)
    }

    pub fn is_paired(&self) -> bool {
        self.0 & 0x01 != 0
    }

    pub fn is_proper_pair(&self) -> bool {
        self.0 & 0x02 != 0
    }

    pub fn is_unmapped(&self) -> bool {
        self.0 & 0x04 != 0
    }

    pub fn is_mate_unmapped(&self) -> bool {
        self.0 & 0x08 != 0
    }

    pub fn is_reverse(&self) -> bool {
        self.0 & 0x10 != 0
    }

    pub fn is_mate_reverse(&self) -> bool {
        self.0 & 0x20 != 0
    }

    pub fn is_read_1(&self) -> bool {
        self.0 & 0x40 != 0
    }

    pub fn is_read_2(&self) -> bool {
        self.0 & 0x80 != 0
    }

    pub fn is_secondary(&self) -> bool {
        self.0 & 0x0100 != 0
    }

    pub fn is_qc_fail(&self) -> bool {
        self.0 & 0x0200 != 0
    }

    pub fn is_dup(&self) -> bool {
        self.0 & 0x0400 != 0
    }

    pub fn is_supplementary(&self) -> bool {
        self.0 & 0x0800 != 0
    }
}

#[derive(Debug)]
pub struct Record {
    ref_id: i32,
    pos: i32,
    mapq: u8,
    pub flag: Flag,
    next_ref_id: i32,
    next_ref_pos: i32,
    pub read_name: String,
    pub cigar: Cigar,
    pub sequence: Sequence,
    pub quality: Quality,
    pub data: Data,
}

impl Record {
    pub fn new(
        ref_id: i32,
        pos: i32,
        mapq: u8,
        flag: Flag,
        next_ref_id: i32,
        next_ref_pos: i32,
        read_name: String,
        cigar: Cigar,
        sequence: Sequence,
        quality: Quality,
        data: Data,
    ) -> Record {
        Record {
            ref_id,
            pos,
            mapq,
            flag,
            next_ref_id,
            next_ref_pos,
            read_name,
            cigar,
            sequence,
            quality,
            data,
        }
    }
}
