mod header;
pub mod records;

use std::{collections::HashMap, io};

use flate2::Compression;
use noodles_fasta as fasta;
use noodles_sam as sam;

pub use self::header::write_header;
use self::records::ExternalDataWriters;
use crate::{
    calculate_normalized_sequence_digest,
    codecs::Encoder,
    container::{
        BlockContentEncoderMap, CompressionHeader, ReferenceSequenceContext,
        block::{self, CompressionMethod, ContentType},
        slice::Header,
    },
    io::{
        BitWriter,
        writer::{Options, Record, container::block::Block},
    },
    record::Flags,
};

pub struct Slice {
    pub header: Header,
    pub core_data_block: Block,
    pub external_data_blocks: Vec<Block>,
}

/// Block content ID used for embedded reference sequences.
///
/// This must not collide with any data series (1-30) or tag encoding content IDs.
const EMBEDDED_REFERENCE_CONTENT_ID: block::ContentId = i32::MAX;

pub(super) fn build_slice(
    reference_sequence_repository: &fasta::Repository,
    options: &Options,
    header: &sam::Header,
    record_counter: u64,
    compression_header: &CompressionHeader,
    records: &mut [Record],
) -> io::Result<Slice> {
    let reference_sequence_context = get_reference_sequence_context(records);

    set_mates(records);

    let (core_data_buf, external_data_bufs) =
        write_records(compression_header, reference_sequence_context, records)?;

    let (core_data_block, mut external_data_blocks) = build_blocks(
        &options.block_content_encoder_map,
        records,
        core_data_buf,
        external_data_bufs,
    )?;

    let mut block_content_ids = vec![core_data_block.content_id];
    block_content_ids.extend(external_data_blocks.iter().map(|block| block.content_id));

    let mut embedded_reference_bases_block_content_id = None;

    let reference_md5 = if options.embed_reference_sequences {
        if let ReferenceSequenceContext::Some(context) = reference_sequence_context {
            let reference_bases =
                get_reference_subsequence(reference_sequence_repository, header, context)?;

            let embedded_ref_block = Block::encode(
                ContentType::ExternalData,
                EMBEDDED_REFERENCE_CONTENT_ID,
                None,
                &reference_bases,
            )?;

            embedded_reference_bases_block_content_id = Some(EMBEDDED_REFERENCE_CONTENT_ID);
            block_content_ids.push(EMBEDDED_REFERENCE_CONTENT_ID);
            external_data_blocks.push(embedded_ref_block);

            // Embedded references don't need MD5 since they ARE the reference.
            None
        } else {
            None
        }
    } else {
        calculate_reference_sequence_md5(
            reference_sequence_repository,
            header,
            reference_sequence_context,
            options.reference_required,
        )?
    };

    let header = Header {
        reference_sequence_context,
        record_count: records.len(),
        record_counter,
        block_count: block_content_ids.len(),
        block_content_ids,
        embedded_reference_bases_block_content_id,
        reference_md5,
        optional_tags: Vec::new(),
    };

    Ok(Slice {
        header,
        core_data_block,
        external_data_blocks,
    })
}

fn get_reference_sequence_context(records: &[Record]) -> ReferenceSequenceContext {
    assert!(!records.is_empty());

    let record = &records[0];
    let mut reference_sequence_context = match (
        record.reference_sequence_id,
        record.alignment_start,
        record.alignment_end(),
    ) {
        (Some(id), Some(start), Some(end)) => ReferenceSequenceContext::some(id, start, end),
        _ => ReferenceSequenceContext::None,
    };

    for record in records.iter().skip(1) {
        reference_sequence_context.update(
            record.reference_sequence_id,
            record.alignment_start,
            record.alignment_end(),
        );
    }

    reference_sequence_context
}

fn set_mates(records: &mut [Record]) {
    assert!(!records.is_empty());

    let mut indices = HashMap::new();
    let mut i = records.len() - 1;

    loop {
        let record = &mut records[i];
        let flags = record.bam_flags;

        if flags.is_segmented() && !flags.is_secondary() {
            let name = record.name.as_ref().map(|name| name.to_owned());

            if let Some(j) = indices.insert(name, i) {
                let mid = i + 1;
                let (left, right) = records.split_at_mut(mid);

                let record = &mut left[i];
                let mate = &mut right[j - mid];

                set_downstream_mate(i, record, j, mate);
            } else {
                set_detached(record);
            }
        } else {
            set_detached(record);
        }

        if i == 0 {
            break;
        }

        i -= 1;
    }
}

fn set_downstream_mate(i: usize, record: &mut Record, j: usize, mate: &mut Record) {
    record.mate_distance = Some(j - i - 1);
    record.cram_flags.insert(Flags::MATE_IS_DOWNSTREAM);
    mate.cram_flags.remove(Flags::IS_DETACHED);
}

fn set_detached(record: &mut Record) {
    record.cram_flags.insert(Flags::IS_DETACHED);
}

#[allow(clippy::type_complexity)]
fn write_records(
    compression_header: &CompressionHeader,
    reference_sequence_context: ReferenceSequenceContext,
    records: &[Record],
) -> io::Result<(Vec<u8>, Vec<(block::ContentId, Vec<u8>)>)> {
    use crate::container::compression_header::data_series_encodings::data_series::STANDARD_DATA_SERIES;

    let mut core_data_writer = BitWriter::default();
    let mut external_data_writers = ExternalDataWriters::default();

    for data_series in STANDARD_DATA_SERIES {
        let block_content_id = block::ContentId::from(*data_series);
        external_data_writers.insert(block_content_id, Vec::new());
    }

    for block_content_id in compression_header.tag_encodings.keys() {
        external_data_writers.insert(*block_content_id, Vec::new());
    }

    let mut writer = records::Writer::new(
        compression_header,
        &mut core_data_writer,
        &mut external_data_writers,
        reference_sequence_context,
    );

    for record in records {
        writer.write_record(record)?;
    }

    Ok((
        core_data_writer.finish()?,
        external_data_writers.into_iter().collect(),
    ))
}

fn build_blocks(
    block_content_encoder_map: &BlockContentEncoderMap,
    records: &[Record],
    core_data_buf: Vec<u8>,
    external_data_bufs: Vec<(block::ContentId, Vec<u8>)>,
) -> io::Result<(Block, Vec<Block>)> {
    use crate::codecs::fqzcomp;

    const CORE_DATA_BLOCK_CONTENT_ID: block::ContentId = 0;
    const DEFAULT_ENCODER: Encoder = Encoder::Gzip(Compression::new(6));

    let core_data_block = Block::encode(
        ContentType::CoreData,
        CORE_DATA_BLOCK_CONTENT_ID,
        block_content_encoder_map.core_data_encoder(),
        &core_data_buf,
    )?;

    let mut all_quality_scores_stored_as_arrays = true;

    for record in records.iter() {
        all_quality_scores_stored_as_arrays = all_quality_scores_stored_as_arrays
            && record.cram_flags.quality_scores_are_stored_as_array();
    }

    let external_data_blocks = external_data_bufs
        .into_iter()
        .filter(|(_, buf)| !buf.is_empty())
        .map(|(block_content_id, buf)| {
            if let Some(encoder) =
                block_content_encoder_map.get_data_series_encoder(block_content_id)
            {
                match encoder {
                    Some(Encoder::Fqzcomp) => {
                        if all_quality_scores_stored_as_arrays {
                            let fqz_records: Vec<_> = records
                                .iter()
                                .map(|r| (r.read_length, r.bam_flags.is_reverse_complemented()))
                                .collect();
                            let data = fqzcomp::encode(&fqz_records, &buf)?;

                            Ok(Block {
                                compression_method: CompressionMethod::Fqzcomp,
                                content_type: ContentType::ExternalData,
                                content_id: block_content_id,
                                uncompressed_size: buf.len(),
                                src: data,
                            })
                        } else {
                            Block::encode(
                                ContentType::ExternalData,
                                block_content_id,
                                Some(&DEFAULT_ENCODER),
                                &buf,
                            )
                        }
                    }
                    _ => Block::encode(ContentType::ExternalData, block_content_id, encoder, &buf),
                }
            } else if let Some(encoder) =
                block_content_encoder_map.get_tag_values_encoders(block_content_id)
            {
                Block::encode(ContentType::ExternalData, block_content_id, encoder, &buf)
            } else {
                Block::encode(
                    ContentType::ExternalData,
                    block_content_id,
                    Some(&DEFAULT_ENCODER),
                    &buf,
                )
            }
        })
        .collect::<io::Result<_>>()?;

    Ok((core_data_block, external_data_blocks))
}

fn calculate_reference_sequence_md5(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    reference_sequence_context: ReferenceSequenceContext,
    reference_required: bool,
) -> io::Result<Option<[u8; 16]>> {
    let ReferenceSequenceContext::Some(context) = reference_sequence_context else {
        return Ok(None);
    };

    let reference_sequence_name = header
        .reference_sequences()
        .get_index(context.reference_sequence_id())
        .map(|(name, _)| name);

    let reference_sequence_name = match reference_sequence_name {
        Some(name) => name,
        None if !reference_required => return Ok(Some([0u8; 16])),
        None => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid reference sequence ID",
            ));
        }
    };

    match reference_sequence_repository.get(reference_sequence_name) {
        Some(result) => {
            let reference_sequence = result?;
            let interval = context.alignment_start()..=context.alignment_end();
            let sequence = &reference_sequence[interval];
            Ok(Some(calculate_normalized_sequence_digest(sequence)))
        }
        None if !reference_required => Ok(Some([0u8; 16])),
        None => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("missing reference sequence: {reference_sequence_name}"),
        )),
    }
}

fn get_reference_subsequence(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    context: crate::container::ReferenceSequenceContextInner,
) -> io::Result<Vec<u8>> {
    let reference_sequence_name = header
        .reference_sequences()
        .get_index(context.reference_sequence_id())
        .map(|(name, _)| name)
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "invalid reference sequence ID")
        })?;

    let reference_sequence = reference_sequence_repository
        .get(reference_sequence_name)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("missing reference sequence: {reference_sequence_name}"),
            )
        })??;

    let interval = context.alignment_start()..=context.alignment_end();
    let sequence = &reference_sequence[interval];

    Ok(sequence.to_vec())
}
