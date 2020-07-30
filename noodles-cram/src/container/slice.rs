mod header;

pub use self::header::Header;

use std::{
    collections::HashMap,
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

    pub fn core_data_block_mut(&mut self) -> &mut Block {
        &mut self.core_data_block
    }

    pub fn add_external_block(&mut self, block: Block) {
        self.external_blocks.push(block);
    }

    pub fn records(&self, compression_header: &CompressionHeader) -> io::Result<Vec<Record>> {
        let core_data = BitReader::new(Cursor::new(self.core_data_block.decompressed_data()));

        let external_data: HashMap<_, _> = self
            .external_blocks
            .iter()
            .map(|block| (block.content_id(), Cursor::new(block.decompressed_data())))
            .collect();

        let mut record_reader = reader::record::Reader::new(
            compression_header,
            core_data,
            external_data,
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
        use std::{
            cell::{Ref, RefCell, RefMut},
            rc::Rc,
        };

        fn set_mate(mut record: RefMut<Record>, mate: Ref<Record>) {
            let mate_bam_flags = mate.bam_bit_flags();

            if mate_bam_flags.is_reverse_complemented() {
                record.bam_bit_flags |=
                    u16::from(sam::record::Flags::MATE_REVERSE_COMPLEMENTED) as i32;
            }

            if mate_bam_flags.is_unmapped() {
                record.bam_bit_flags |= u16::from(sam::record::Flags::MATE_UNMAPPED) as i32;
            }

            record.next_fragment_reference_sequence_id = mate.reference_id;
            record.next_mate_alignment_start = mate.alignment_start;
        }

        let records: Vec<_> = records
            .into_iter()
            .map(|r| Rc::new(RefCell::new(r)))
            .collect();

        let mut mate_indicies = vec![None; records.len()];

        for (i, record_cell) in records.iter().enumerate() {
            let record = record_cell.borrow();
            let cram_flags = record.cram_bit_flags();

            if cram_flags.has_mate_downstream() {
                let distance_to_next_fragment = record.distance_to_next_fragment as usize;
                let mate_index = i + distance_to_next_fragment + 1;
                mate_indicies[i] = Some(mate_index);
            }
        }

        for (i, record_cell) in records.iter().enumerate() {
            if mate_indicies[i].is_none() {
                continue;
            }

            let mut record = record_cell.borrow_mut();
            let mut j = i;

            while let Some(mate_index) = mate_indicies[j] {
                let mate = records[mate_index].borrow();
                set_mate(record, mate);
                record = records[mate_index].borrow_mut();
                j = mate_index;
            }

            let mate = record_cell.borrow();
            set_mate(record, mate);

            // TODO: calculate template size
        }

        records
            .into_iter()
            .map(|r| Rc::try_unwrap(r).unwrap().into_inner())
            .collect()
    }
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
