mod builder;

pub use builder::Builder;

use std::io;

use noodles_sam as sam;

use crate::container::ReferenceSequenceId;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Header {
    reference_sequence_id: ReferenceSequenceId,
    alignment_start: Option<sam::record::Position>,
    alignment_span: usize,
    record_count: usize,
    record_counter: i64,
    block_count: usize,
    block_content_ids: Vec<i32>,
    embedded_reference_bases_block_content_id: Option<i32>,
    reference_md5: [u8; 16],
    optional_tags: Vec<u8>,
}

impl Header {
    pub fn builder() -> Builder {
        Builder::default()
    }

    pub fn reference_sequence_id(&self) -> ReferenceSequenceId {
        self.reference_sequence_id
    }

    pub fn alignment_start(&self) -> Option<sam::record::Position> {
        self.alignment_start
    }

    pub fn alignment_span(&self) -> usize {
        self.alignment_span
    }

    pub fn alignment_end(&self) -> Option<io::Result<sam::record::Position>> {
        self.alignment_start().map(|start| {
            let start = i32::from(start);
            let len = self.alignment_span() as i32;
            let end = start + len + 1;

            sam::record::Position::try_from(end)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }

    pub fn record_count(&self) -> usize {
        self.record_count
    }

    pub fn record_counter(&self) -> i64 {
        self.record_counter
    }

    pub fn block_count(&self) -> usize {
        self.block_count
    }

    pub fn block_content_ids(&self) -> &[i32] {
        &self.block_content_ids
    }

    pub fn embedded_reference_bases_block_content_id(&self) -> Option<i32> {
        self.embedded_reference_bases_block_content_id
    }

    pub fn reference_md5(&self) -> &[u8] {
        &self.reference_md5
    }

    pub fn optional_tags(&self) -> &[u8] {
        &self.optional_tags
    }
}
