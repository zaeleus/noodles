use crate::{
    num::{Itf8, Ltf8},
    CompressionHeader, Slice,
};

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
}
