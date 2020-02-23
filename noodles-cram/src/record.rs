#[derive(Clone, Copy, Debug)]
pub struct CramBitFlags(i32);

impl CramBitFlags {
    pub fn are_quality_scores_stored_as_array(self) -> bool {
        self.0 & 0x01 != 0
    }

    pub fn is_detached(self) -> bool {
        self.0 & 0x02 != 0
    }

    pub fn has_mate_downstream(self) -> bool {
        self.0 & 0x04 != 0
    }

    pub fn decode_sequence_as_unknown(self) -> bool {
        self.0 & 0x08 != 0
    }
}

#[derive(Clone, Debug, Default)]
pub struct Record {
    pub bam_bit_flags: i32,
    pub cram_bit_flags: i32,
    pub reference_id: i32,
    pub read_length: i32,
    pub alignment_start: i32,
    pub read_group: i32,
    pub read_name: Vec<u8>,
    pub next_mate_bit_flags: i32,
    pub next_fragment_reference_sequence_id: i32,
    pub next_mate_alignment_start: i32,
    pub template_size: i32,
    pub distance_to_next_fragment: i32,
}

impl Record {
    pub fn cram_bit_flags(&self) -> CramBitFlags {
        CramBitFlags(self.cram_bit_flags)
    }
}
