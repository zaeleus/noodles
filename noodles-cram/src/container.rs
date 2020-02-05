use crate::{
    num::{Itf8, Ltf8},
    CompressionHeader,
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
        }
    }

    pub fn len(&self) -> i32 {
        self.length
    }
}

#[derive(Debug, Default)]
pub struct Container {
    header: Header,
    compression_header: CompressionHeader,
    blocks: Vec<u8>,
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

    pub fn blocks(&self) -> &[u8] {
        &self.blocks
    }

    pub fn blocks_mut(&mut self) -> &mut Vec<u8> {
        &mut self.blocks
    }
}
