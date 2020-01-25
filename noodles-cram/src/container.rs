use crate::num::{Itf8, Ltf8};

#[derive(Debug, Default)]
pub struct Header {
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
}

#[derive(Debug, Default)]
pub struct Container {
    header: Header,
    blocks: Vec<u8>,
}

impl Container {
    pub fn new(header: Header, blocks: Vec<u8>) -> Self {
        Self { header, blocks }
    }

    pub fn header(&self) -> &Header {
        &self.header
    }
}
