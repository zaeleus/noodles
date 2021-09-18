pub(crate) mod builder;
pub(crate) mod header;

pub use self::{builder::Builder, header::Header};

use std::{
    collections::HashMap,
    convert::TryFrom,
    io::{self, Cursor},
};

use noodles_sam as sam;

use super::CompressionHeader;
use crate::{container::Block, BitReader, Record};

/// A CRAM data container slice.
///
/// A slice contains a header, a core data block, and one or more external blocks. This is where
/// the CRAM records are stored.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Slice {
    header: Header,
    core_data_block: Block,
    external_blocks: Vec<Block>,
}

impl Slice {
    pub(crate) fn builder() -> Builder {
        Builder::default()
    }

    pub(crate) fn new(header: Header, core_data_block: Block, external_blocks: Vec<Block>) -> Self {
        Self {
            header,
            core_data_block,
            external_blocks,
        }
    }

    pub(crate) fn header(&self) -> &Header {
        &self.header
    }

    pub(crate) fn core_data_block(&self) -> &Block {
        &self.core_data_block
    }

    pub(crate) fn external_blocks(&self) -> &[Block] {
        &self.external_blocks
    }

    pub fn records(&self, compression_header: &CompressionHeader) -> io::Result<Vec<Record>> {
        let core_data_reader = self
            .core_data_block
            .decompressed_data()
            .map(Cursor::new)
            .map(BitReader::new)?;

        let mut external_data_readers = HashMap::with_capacity(self.external_blocks().len());

        for block in self.external_blocks() {
            let reader = block.decompressed_data().map(Cursor::new)?;
            external_data_readers.insert(block.content_id(), reader);
        }

        let mut record_reader = crate::reader::record::Reader::new(
            compression_header,
            core_data_reader,
            external_data_readers,
            self.header.reference_sequence_id(),
            self.header.alignment_start(),
        );

        let record_count = self.header().record_count();
        let mut records = Vec::with_capacity(record_count);

        let start_id = self.header().record_counter();
        let end_id = start_id + (record_count as i64);

        for id in start_id..end_id {
            let mut record = record_reader.read_record()?;
            record.id = id;
            records.push(record);
        }

        Ok(records)
    }

    pub fn resolve_mates(&self, records: Vec<Record>) -> Vec<Record> {
        use std::cell::RefCell;

        let mut mate_indices = vec![None; records.len()];

        for (i, record) in records.iter().enumerate() {
            let flags = record.flags();

            if flags.has_mate_downstream() {
                let distance_to_next_fragment = record.distance_to_next_fragment() as usize;
                let mate_index = i + distance_to_next_fragment + 1;
                mate_indices[i] = Some(mate_index);
            }
        }

        let records: Vec<_> = records.into_iter().map(RefCell::new).collect();

        for (i, record_cell) in records.iter().enumerate() {
            if mate_indices[i].is_none() {
                continue;
            }

            let mut record = record_cell.borrow_mut();

            if record.read_name.is_empty() {
                let read_name = record.id().to_string().into_bytes();
                record.read_name.extend(read_name);
            }

            let mut j = i;

            while let Some(mate_index) = mate_indices[j] {
                let mut mate = records[mate_index].borrow_mut();
                set_mate(&mut record, &mut mate);
                record = mate;
                j = mate_index;
            }

            let mut mate = record_cell.borrow_mut();
            set_mate(&mut record, &mut mate);

            let template_size = calculate_template_size(&record, &mate);
            record.template_size = template_size;
            mate.template_size = -template_size;
        }

        records.into_iter().map(|r| r.into_inner()).collect()
    }
}

fn set_mate(mut record: &mut Record, mate: &mut Record) {
    let mate_bam_flags = mate.bam_flags();

    if mate_bam_flags.is_reverse_complemented() {
        record.bam_bit_flags |= sam::record::Flags::MATE_REVERSE_COMPLEMENTED;
    }

    if mate_bam_flags.is_unmapped() {
        record.bam_bit_flags |= sam::record::Flags::MATE_UNMAPPED;
    }

    if mate.read_name().is_empty() {
        mate.read_name.extend(record.read_name.iter());
    }

    record.next_fragment_reference_sequence_id = mate.reference_sequence_id();
    record.next_mate_alignment_start = mate.alignment_start();
}

fn calculate_template_size(record: &Record, mate: &Record) -> i32 {
    let start = record.alignment_start().map(i32::from).unwrap_or_default();
    let end = mate.alignment_end();
    end - start + 1
}

impl TryFrom<&[Block]> for Slice {
    type Error = io::Error;

    fn try_from(blocks: &[Block]) -> Result<Self, Self::Error> {
        let mut it = blocks.iter();

        let header_block = it
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing header block"))?;
        let data = header_block.decompressed_data()?;
        let mut reader = &data[..];
        let header = crate::reader::data_container::slice::read_header(&mut reader)?;

        let core_data_block = it
            .next()
            .cloned()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing core data block"))?;

        let external_block_count = (header.block_count() - 1) as usize;
        let external_blocks: Vec<_> = it.take(external_block_count).cloned().collect();

        assert_eq!(external_block_count, external_blocks.len());

        Ok(Slice::new(header, core_data_block, external_blocks))
    }
}
