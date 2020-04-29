use crate::{
    num::{Itf8, Ltf8},
    CompressionHeader, Slice,
};

const EOF_LEN: i32 = 15;
const EOF_REFERENCE_SEQUENCE_ID: Itf8 = -1;
const EOF_START_POSITION: Itf8 = 4_542_278;
const EOF_N_BLOCKS: Itf8 = 1;
const EOF_CRC32: u32 = 0x4f_d9_bd_05;

#[derive(Debug, Default)]
pub struct Header {
    length: i32,
    reference_sequence_id: Itf8,
    start_position: Itf8,
    alignment_span: Itf8,
    n_records: Itf8,
    record_counter: Ltf8,
    bases: Ltf8,
    n_blocks: Itf8,
    landmarks: Vec<Itf8>,
    crc32: u32,
}

impl Header {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        length: i32,
        reference_sequence_id: Itf8,
        start_position: Itf8,
        alignment_span: Itf8,
        n_records: Itf8,
        record_counter: Ltf8,
        bases: Ltf8,
        n_blocks: Itf8,
        landmarks: Vec<Itf8>,
        crc32: u32,
    ) -> Self {
        Self {
            length,
            reference_sequence_id,
            start_position,
            alignment_span,
            n_records,
            record_counter,
            bases,
            n_blocks,
            landmarks,
            crc32,
        }
    }

    pub fn len(&self) -> i32 {
        self.length
    }

    pub fn reference_sequence_id(&self) -> Itf8 {
        self.reference_sequence_id
    }

    pub fn is_eof(&self) -> bool {
        self.length == EOF_LEN
            && self.reference_sequence_id == EOF_REFERENCE_SEQUENCE_ID
            && self.start_position == EOF_START_POSITION
            && self.alignment_span == 0
            && self.n_records == 0
            && self.record_counter == 0
            && self.bases == 0
            && self.n_blocks == EOF_N_BLOCKS
            && self.landmarks.is_empty()
            && self.crc32 == EOF_CRC32
    }

    pub fn landmarks(&self) -> &[Itf8] {
        &self.landmarks
    }
}

#[derive(Debug, Default)]
pub struct Container {
    header: Header,
    compression_header: CompressionHeader,
    slices: Vec<Slice>,
}

impl Container {
    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn header_mut(&mut self) -> &mut Header {
        &mut self.header
    }

    pub fn compression_header(&self) -> &CompressionHeader {
        &self.compression_header
    }

    pub fn compression_header_mut(&mut self) -> &mut CompressionHeader {
        &mut self.compression_header
    }

    pub fn slices(&self) -> &[Slice] {
        &self.slices
    }

    pub fn slices_mut(&mut self) -> &mut Vec<Slice> {
        &mut self.slices
    }

    pub fn add_slice(&mut self, slice: Slice) {
        self.slices.push(slice);
    }

    pub fn clear(&mut self) {
        self.header = Default::default();
        self.compression_header = Default::default();
        self.slices.clear();
    }

    pub fn is_eof(&self) -> bool {
        self.header.is_eof()
    }
}
