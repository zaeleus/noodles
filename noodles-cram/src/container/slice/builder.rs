use std::{cmp, collections::HashMap, io};

use md5::{Digest, Md5};
use noodles_fasta as fasta;

use crate::{
    container::{
        block::{self, CompressionMethod},
        Block, ReferenceSequenceId,
    },
    data_container::{compression_header::data_series_encoding_map::DataSeries, CompressionHeader},
    writer, BitWriter, Record,
};

use super::{Header, Slice};

use noodles_bam as bam;

const CORE_DATA_BLOCK_CONTENT_ID: i32 = 0;
const MAX_RECORD_COUNT: usize = 2560;

#[derive(Debug, Default)]
pub struct Builder {
    records: Vec<Record>,
    reference_sequence_id: Option<Option<bam::record::ReferenceSequenceId>>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AddRecordError {
    ReferenceSequenceIdMismatch(Record),
    SliceFull(Record),
}

impl Builder {
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn add_record(&mut self, record: Record) -> Result<&Record, AddRecordError> {
        if self.records.len() >= MAX_RECORD_COUNT {
            return Err(AddRecordError::SliceFull(record));
        }

        let record_reference_sequence_id = record.reference_sequence_id();

        if self.reference_sequence_id.is_none() {
            self.reference_sequence_id = Some(record_reference_sequence_id);
        }

        match self.reference_sequence_id.unwrap() {
            Some(slice_reference_sequence_id) => {
                match record_reference_sequence_id.map(i32::from) {
                    Some(id) => {
                        if i32::from(slice_reference_sequence_id) == id {
                            self.records.push(record);
                            Ok(self.records.last().unwrap())
                        } else {
                            Err(AddRecordError::ReferenceSequenceIdMismatch(record))
                        }
                    }
                    None => Err(AddRecordError::ReferenceSequenceIdMismatch(record)),
                }
            }
            None => match record_reference_sequence_id {
                Some(_) => Err(AddRecordError::ReferenceSequenceIdMismatch(record)),
                None => {
                    self.records.push(record);
                    Ok(self.records.last().unwrap())
                }
            },
        }
    }

    pub fn build(
        self,
        reference_sequences: &[fasta::Record],
        compression_header: &CompressionHeader,
        record_counter: i64,
    ) -> io::Result<Slice> {
        let reference_sequence_id = match self.reference_sequence_id.unwrap() {
            Some(id) => ReferenceSequenceId::Some(i32::from(id)),
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

        let mut slice_alignment_start = i32::MAX;
        let mut slice_alignment_end = 1;

        for record in &self.records {
            slice_alignment_start = cmp::min(slice_alignment_start, record.alignment_start());
            slice_alignment_end = cmp::max(slice_alignment_end, record.alignment_end());

            record_writer.write_record(record)?;
        }

        let core_data_block = core_data_writer.finish().and_then(|buf| {
            Block::builder()
                .set_content_type(block::ContentType::CoreData)
                .set_content_id(CORE_DATA_BLOCK_CONTENT_ID)
                .compress_and_set_data(buf, CompressionMethod::Gzip)
                .map(|builder| builder.build())
        })?;

        let mut block_content_ids = vec![CORE_DATA_BLOCK_CONTENT_ID];

        let external_blocks: Vec<_> = external_data_writers
            .into_iter()
            .filter(|(_, buf)| !buf.is_empty())
            .map(|(block_content_id, buf)| {
                Block::builder()
                    .set_content_type(block::ContentType::ExternalData)
                    .set_content_id(block_content_id)
                    .compress_and_set_data(buf, CompressionMethod::Gzip)
                    .map(|builder| builder.build())
            })
            .collect::<Result<_, _>>()?;

        for block in &external_blocks {
            block_content_ids.push(block.content_id());
        }

        let reference_md5 = if let ReferenceSequenceId::Some(id) = reference_sequence_id {
            let reference_sequence = reference_sequences
                .get(id as usize)
                .map(|record| record.sequence())
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "missing reference sequence")
                })?;

            let start = (slice_alignment_start - 1) as usize;
            let end = (slice_alignment_end - 1) as usize;

            let mut hasher = Md5::new();
            hasher.update(&reference_sequence[start..=end]);
            <[u8; 16]>::from(hasher.finalize())
        } else {
            [0; 16]
        };

        let slice_alignment_span = slice_alignment_end - slice_alignment_start + 1;

        let header = Header::builder()
            .set_reference_sequence_id(reference_sequence_id)
            .set_alignment_start(slice_alignment_start)
            .set_alignment_span(slice_alignment_span)
            .set_record_count(self.records.len() as i32)
            .set_record_counter(record_counter)
            // external blocks + core data block
            .set_block_count(external_blocks.len() + 1)
            .set_block_content_ids(block_content_ids)
            .set_reference_md5(reference_md5)
            .build();

        Ok(Slice::new(header, core_data_block, external_blocks))
    }
}
