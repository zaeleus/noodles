use std::io::{self, Write};

use super::{fmt_display_field, write_other_fields};
use crate::header::record::value::{
    map::{read_group::tag, ReadGroup},
    Map,
};

pub(crate) fn write_read_group<W>(
    writer: &mut W,
    id: &[u8],
    read_group: &Map<ReadGroup>,
) -> io::Result<()>
where
    W: Write,
{
    use crate::io::writer::header::record::{value::map::write_separator, write_delimiter};

    write_delimiter(writer)?;
    writer.write_all(tag::ID.as_ref())?;
    write_separator(writer)?;
    writer.write_all(id)?;

    if let Some(barcode) = read_group.barcode() {
        fmt_display_field(writer, tag::BARCODE, barcode)?;
    }

    if let Some(sequencing_center) = read_group.sequencing_center() {
        fmt_display_field(writer, tag::SEQUENCING_CENTER, sequencing_center)?;
    }

    if let Some(description) = read_group.description() {
        fmt_display_field(writer, tag::DESCRIPTION, description)?;
    }

    if let Some(produced_at) = read_group.produced_at() {
        fmt_display_field(writer, tag::PRODUCED_AT, produced_at)?;
    }

    if let Some(flow_order) = read_group.flow_order() {
        fmt_display_field(writer, tag::FLOW_ORDER, flow_order)?;
    }

    if let Some(key_sequence) = read_group.key_sequence() {
        fmt_display_field(writer, tag::KEY_SEQUENCE, key_sequence)?;
    }

    if let Some(library) = read_group.library() {
        fmt_display_field(writer, tag::LIBRARY, library)?;
    }

    if let Some(program) = read_group.program() {
        fmt_display_field(writer, tag::PROGRAM, program)?;
    }

    if let Some(predicted_median_insert_size) = read_group.predicted_median_insert_size() {
        fmt_display_field(
            writer,
            tag::PREDICTED_MEDIAN_INSERT_SIZE,
            predicted_median_insert_size,
        )?;
    }

    if let Some(platform) = read_group.platform() {
        fmt_display_field(writer, tag::PLATFORM, platform)?;
    }

    if let Some(platform_model) = read_group.platform_model() {
        fmt_display_field(writer, tag::PLATFORM_MODEL, platform_model)?;
    }

    if let Some(platform_unit) = read_group.platform_unit() {
        fmt_display_field(writer, tag::PLATFORM_UNIT, platform_unit)?;
    }

    if let Some(sample) = read_group.sample() {
        fmt_display_field(writer, tag::SAMPLE, sample)?;
    }

    write_other_fields(writer, read_group.other_fields())?;

    Ok(())
}
