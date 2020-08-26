pub mod builder;
pub mod header;

pub use self::{builder::Builder, header::Header};

use std::{
    convert::TryFrom,
    io::{self, Cursor},
};

use noodles_sam as sam;

use crate::{reader, BitReader, Record};

use super::{Block, CompressionHeader};

#[derive(Debug)]
pub struct Slice {
    header: Header,
    core_data_block: Block,
    external_blocks: Vec<Block>,
}

impl Slice {
    pub fn builder() -> Builder {
        Builder::default()
    }

    pub fn new(header: Header, core_data_block: Block, external_blocks: Vec<Block>) -> Self {
        Self {
            header,
            core_data_block,
            external_blocks,
        }
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn core_data_block(&self) -> &Block {
        &self.core_data_block
    }

    pub fn external_blocks(&self) -> &[Block] {
        &self.external_blocks
    }

    pub fn records(&self, compression_header: &CompressionHeader) -> io::Result<Vec<Record>> {
        let core_data_reader =
            BitReader::new(Cursor::new(self.core_data_block.decompressed_data()));

        let external_data_readers = self
            .external_blocks
            .iter()
            .map(|block| (block.content_id(), Cursor::new(block.decompressed_data())))
            .collect();

        let mut record_reader = reader::record::Reader::new(
            compression_header,
            core_data_reader,
            external_data_readers,
            self.header.reference_sequence_id(),
            self.header.alignment_start(),
        );

        let record_counter = self.header().record_counter();
        let records_len = self.header().record_count() as usize;
        let mut records = Vec::with_capacity(records_len);

        for i in 0..records_len {
            let mut record = Record::default();
            record.id = record_counter + (i as i64);
            record_reader.read_record(&mut record)?;
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
                let distance_to_next_fragment = record.distance_to_next_fragment as usize;
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
                let read_name = record.id.to_string().into_bytes();
                record.read_name.extend(read_name);
            }

            let mut j = i;

            while let Some(mate_index) = mate_indices[j] {
                let mut mate = records[mate_index].borrow_mut();
                set_mate(&mut record, &mut mate);
                record = records[mate_index].borrow_mut();
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

    if mate.read_name.is_empty() {
        mate.read_name.extend(record.read_name.iter());
    }

    record.next_fragment_reference_sequence_id = mate.reference_sequence_id();
    record.next_mate_alignment_start = mate.alignment_start;
}

fn calculate_template_size(record: &Record, mate: &Record) -> i32 {
    let start = record.alignment_start();
    let end = mate.alignment_end();
    end - start + 1
}

impl TryFrom<&[Block]> for Slice {
    type Error = io::Error;

    fn try_from(blocks: &[Block]) -> Result<Self, Self::Error> {
        let data = blocks[0].decompressed_data();
        let mut reader = &data[..];
        let header = reader::slice::read_header(&mut reader)?;

        let core_data_block = blocks[1].clone();

        let external_blocks_len = header.block_count() as usize;
        let mut external_blocks = Vec::with_capacity(external_blocks_len);

        for i in 0..external_blocks_len {
            external_blocks.push(blocks[i + 1].clone());
        }

        Ok(Slice::new(header, core_data_block, external_blocks))
    }
}
