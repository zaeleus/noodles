use std::{cmp, collections::HashMap, io};

use md5::{Digest, Md5};
use noodles_fasta as fasta;
use noodles_sam as sam;

use crate::{
    container::{
        block::{self, CompressionMethod},
        Block, ReferenceSequenceId,
    },
    data_container::{
        compression_header::{data_series_encoding_map::DataSeries, SubstitutionMatrix},
        CompressionHeader,
    },
    record::{Feature, Features},
    writer, BitWriter, Record,
};

use super::{Header, Slice};

use noodles_bam as bam;

const CORE_DATA_BLOCK_CONTENT_ID: i32 = 0;
const MAX_RECORD_COUNT: usize = 10240;

#[derive(Debug, Default)]
pub struct Builder {
    records: Vec<Record>,
    slice_reference_sequence_id: Option<bam::record::ReferenceSequenceId>,
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

        if self.records.is_empty() {
            self.slice_reference_sequence_id = record.reference_sequence_id();
        }

        if record.reference_sequence_id() == self.slice_reference_sequence_id {
            self.records.push(record);
            Ok(self.records.last().unwrap())
        } else {
            Err(AddRecordError::ReferenceSequenceIdMismatch(record))
        }
    }

    pub fn build(
        mut self,
        reference_sequences: &[fasta::Record],
        compression_header: &CompressionHeader,
        record_counter: i64,
    ) -> io::Result<Slice> {
        let reference_sequence_id = match self.slice_reference_sequence_id {
            Some(id) => {
                // FIXME
                ReferenceSequenceId::Some(usize::from(id) as i32)
            }
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
        let mut slice_alignment_end = 0;

        for record in &mut self.records {
            if let ReferenceSequenceId::Some(id) = reference_sequence_id {
                if let Some(alignment_start) = record.alignment_start() {
                    let reference_sequence = reference_sequences
                        .get(id as usize)
                        .map(|record| record.sequence())
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidInput,
                                "missing reference sequence",
                            )
                        })?;

                    update_substitution_codes(
                        reference_sequence.as_ref(),
                        compression_header.preservation_map().substitution_matrix(),
                        alignment_start,
                        &record.bases,
                        &mut record.features,
                    );
                }
            }

            let record_alignment_start =
                record.alignment_start().map(i32::from).unwrap_or_default();

            slice_alignment_start = cmp::min(slice_alignment_start, record_alignment_start);

            let record_alignment_end = record
                .alignment_end()
                .transpose()?
                .map(i32::from)
                .unwrap_or_default();

            slice_alignment_end = cmp::max(slice_alignment_end, record_alignment_end);

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
            hasher.update(&reference_sequence.as_ref()[start..=end]);
            <[u8; 16]>::from(hasher.finalize())
        } else {
            [0; 16]
        };

        let mut builder = Header::builder()
            .set_reference_sequence_id(reference_sequence_id)
            .set_record_count(self.records.len())
            .set_record_counter(record_counter)
            // external blocks + core data block
            .set_block_count(external_blocks.len() + 1)
            .set_block_content_ids(block_content_ids)
            .set_reference_md5(reference_md5);

        if reference_sequence_id.is_some() {
            let slice_alignment_span = slice_alignment_end - slice_alignment_start + 1;

            let slice_alignment_start = sam::record::Position::try_from(slice_alignment_start)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

            builder = builder
                .set_alignment_start(slice_alignment_start)
                .set_alignment_span(slice_alignment_span);
        }

        let header = builder.build();

        Ok(Slice::new(header, core_data_block, external_blocks))
    }
}

fn update_substitution_codes(
    reference_sequence: &[u8],
    substitution_matrix: &SubstitutionMatrix,
    alignment_start: sam::record::Position,
    read_bases: &[u8],
    features: &mut Features,
) {
    use crate::data_container::compression_header::preservation_map::substitution_matrix::Base;

    let mut codes = Vec::new();

    for ((mut reference_position, mut read_position), feature) in
        features.with_positions(alignment_start)
    {
        reference_position -= 1;
        read_position -= 1;

        if let Feature::Substitution(..) = feature {
            let base = char::from(reference_sequence[reference_position]);
            let reference_base = Base::try_from(base).unwrap_or_default();

            let base = char::from(read_bases[read_position]);
            let read_base = Base::try_from(base).unwrap_or_default();

            let code = substitution_matrix.find_code(reference_base, read_base);
            codes.push(code);
        }
    }

    let mut i = 0;

    for feature in features.iter_mut() {
        if let Feature::Substitution(pos, _) = feature {
            *feature = Feature::Substitution(*pos, codes[i]);
            i += 1;
        }
    }
}
