mod substitution_matrix;
mod tag_sets;

use std::io::{self, Write};

use self::{
    substitution_matrix::{build_substitution_matrix, write_substitution_matrix},
    tag_sets::{build_tag_sets, write_tag_sets},
};
use crate::{
    container::compression_header::{PreservationMap, preservation_map::Key},
    file_definition::Version,
    io::writer::{
        Options, Record,
        collections::write_array,
        num::{write_int, write_u8},
    },
};

pub(super) fn write_preservation_map<W>(
    writer: &mut W,
    preservation_map: &PreservationMap,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    let buf = encode(preservation_map, version)?;
    write_array(writer, version, &buf)
}

fn encode(preservation_map: &PreservationMap, version: Version) -> io::Result<Vec<u8>> {
    let mut buf = Vec::new();
    encode_inner(&mut buf, preservation_map, version)?;
    Ok(buf)
}

fn encode_inner<W>(
    writer: &mut W,
    preservation_map: &PreservationMap,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    let map_length: i32 = if version >= Version::V4_0 { 6 } else { 5 };

    write_int(writer, version, map_length)?;

    write_key(writer, Key::RecordsHaveNames)?;
    write_bool(writer, preservation_map.records_have_names())?;

    write_key(writer, Key::AlignmentStartsAreDeltas)?;
    write_bool(writer, preservation_map.alignment_starts_are_deltas())?;

    write_key(writer, Key::ExternalReferenceSequenceIsRequired)?;
    write_bool(
        writer,
        preservation_map.external_reference_sequence_is_required(),
    )?;

    write_key(writer, Key::SubstitutionMatrix)?;
    write_substitution_matrix(writer, preservation_map.substitution_matrix())?;

    write_key(writer, Key::TagSets)?;
    write_tag_sets(writer, preservation_map.tag_sets(), version)?;

    if version >= Version::V4_0 {
        write_key(writer, Key::QualityScoreOrientation)?;
        write_quality_score_orientation(writer, preservation_map.qs_seq_orient())?;
    }

    Ok(())
}

fn write_quality_score_orientation<W>(writer: &mut W, qs_seq_orient: bool) -> io::Result<()>
where
    W: Write,
{
    // CRAM 4.0: 0 = original/sequencing orientation, 1 = alignment orientation.
    let value: u8 = if qs_seq_orient { 1 } else { 0 };
    write_u8(writer, value)
}

fn write_key<W>(writer: &mut W, key: Key) -> io::Result<()>
where
    W: Write,
{
    let data = <[u8; 2]>::from(key);
    writer.write_all(&data)
}

fn write_bool<W>(writer: &mut W, value: bool) -> io::Result<()>
where
    W: Write,
{
    const FALSE: u8 = 0x00;
    const TRUE: u8 = 0x01;

    if value {
        write_u8(writer, TRUE)
    } else {
        write_u8(writer, FALSE)
    }
}

pub(super) fn build_preservation_map(options: &Options, records: &[Record]) -> PreservationMap {
    // § 8.4 Compression header block (2020-06-22): "The boolean values are optional, defaulting to
    // true when absent, although it is recommended to explicitly set them."
    PreservationMap {
        records_have_names: options.preserve_read_names,
        alignment_starts_are_deltas: options.encode_alignment_start_positions_as_deltas,
        external_reference_sequence_is_required: options.reference_required
            && !options.embed_reference_sequences,
        substitution_matrix: build_substitution_matrix(records),
        tag_sets: build_tag_sets(records),
        qs_seq_orient: options.qs_seq_orient,
    }
}
