#[derive(Clone, Debug, Default)]
pub struct Record {
    pub bam_bit_flags: i32,
    pub cram_bit_flags: i32,
    pub reference_id: i32,
    pub read_length: i32,
    pub alignment_start: i32,
    pub read_group: i32,
    pub read_name: Vec<u8>,
}
