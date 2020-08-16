use std::{cmp, collections::HashMap, io};

use crate::{
    container::{
        block, compression_header::data_series_encoding_map::DataSeries, Block, CompressionHeader,
        ReferenceSequenceId,
    },
    writer, BitWriter, Record,
};

use super::{header::EmbeddedReferenceBasesBlockContentId, Header, Slice};

use noodles_bam as bam;

const CORE_DATA_BLOCK_CONTENT_ID: i32 = 0;

#[derive(Debug, Default)]
pub struct Builder {
    records: Vec<Record>,
    reference_sequence_id: Option<bam::record::ReferenceSequenceId>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AddRecordError {
    ReferenceSequenceIdMismatch(Record),
}

impl Builder {
    pub fn add_record(&mut self, record: Record) -> Result<&Record, AddRecordError> {
        let record_reference_sequence_id =
            bam::record::ReferenceSequenceId::from(record.reference_id);

        if self.reference_sequence_id.is_none() {
            self.reference_sequence_id = Some(record_reference_sequence_id);
        }

        match *self.reference_sequence_id.unwrap() {
            Some(slice_reference_sequence_id) => match *record_reference_sequence_id {
                Some(id) => {
                    if slice_reference_sequence_id == id {
                        self.records.push(record);
                        Ok(self.records.last().unwrap())
                    } else {
                        Err(AddRecordError::ReferenceSequenceIdMismatch(record))
                    }
                }
                None => Err(AddRecordError::ReferenceSequenceIdMismatch(record)),
            },
            None => match *record_reference_sequence_id {
                Some(_) => Err(AddRecordError::ReferenceSequenceIdMismatch(record)),
                None => {
                    self.records.push(record);
                    Ok(self.records.last().unwrap())
                }
            },
        }
    }

    pub fn build(self, compression_header: &CompressionHeader) -> io::Result<Slice> {
        let reference_sequence_id = match *self.reference_sequence_id.unwrap() {
            Some(id) => ReferenceSequenceId::Some(id),
            None => ReferenceSequenceId::None,
        };

        let alignment_start = self
            .records
            .first()
            .map(|r| r.alignment_start())
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "no records in builder"))?;

        let mut core_data_writer = BitWriter::new(Vec::new());

        let mut external_data_writers = HashMap::new();

        for i in 0..DataSeries::LEN {
            let block_content_id = (i + 1) as i32;
            external_data_writers.insert(block_content_id, Vec::new());
        }

        for &block_content_id in compression_header.tag_encoding_map().keys() {
            external_data_writers.insert(block_content_id, Vec::new());
        }

        let mut record_writer = writer::record::Writer::new(
            compression_header,
            &mut core_data_writer,
            &mut external_data_writers,
            reference_sequence_id,
            alignment_start,
        );

        let mut min_alignment_start = i32::MAX;
        let mut max_alignment_end = 1;

        for record in &self.records {
            min_alignment_start = cmp::min(min_alignment_start, record.alignment_start());
            max_alignment_end = cmp::max(max_alignment_end, record.alignment_end());

            record_writer.write_record(record)?;
        }

        let core_data_block = core_data_writer.finish().map(|buf| {
            Block::new(
                block::CompressionMethod::None,
                block::ContentType::CoreData,
                CORE_DATA_BLOCK_CONTENT_ID,
                buf.len() as i32,
                buf,
                0,
            )
        })?;

        let mut block_content_ids = vec![CORE_DATA_BLOCK_CONTENT_ID];

        for &block_content_id in external_data_writers.keys() {
            block_content_ids.push(block_content_id);
        }

        let external_blocks: Vec<_> = external_data_writers
            .into_iter()
            .filter(|(_, buf)| !buf.is_empty())
            .map(|(block_content_id, buf)| {
                Block::new(
                    block::CompressionMethod::None,
                    block::ContentType::ExternalData,
                    block_content_id,
                    buf.len() as i32,
                    buf,
                    0,
                )
            })
            .collect();

        // TODO
        let header = Header::new(
            reference_sequence_id,
            min_alignment_start,
            max_alignment_end,
            self.records.len() as i32,
            0,
            (external_blocks.len() + 1) as i32,
            block_content_ids,
            EmbeddedReferenceBasesBlockContentId::default(),
            Default::default(),
            Vec::new(),
        );

        Ok(Slice::new(header, core_data_block, external_blocks))
    }
}
