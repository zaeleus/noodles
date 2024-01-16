mod length;
mod name;

use std::io::{self, Write};

use self::{length::write_length_field, name::write_name_field};
use super::{write_field, write_other_fields};
use crate::header::record::value::{
    map::{reference_sequence::tag, ReferenceSequence},
    Map,
};

pub(crate) fn write_reference_sequence<W>(
    writer: &mut W,
    name: &[u8],
    reference_sequence: &Map<ReferenceSequence>,
) -> io::Result<()>
where
    W: Write,
{
    write_name_field(writer, name)?;
    write_length_field(writer, reference_sequence.length())?;

    if let Some(alternative_locus) = reference_sequence.alternative_locus() {
        write_field(writer, tag::ALTERNATIVE_LOCUS, alternative_locus)?;
    }

    if let Some(alternative_names) = reference_sequence.alternative_names() {
        write_field(writer, tag::ALTERNATIVE_NAMES, alternative_names)?;
    }

    if let Some(assembly_id) = reference_sequence.assembly_id() {
        write_field(writer, tag::ASSEMBLY_ID, assembly_id)?;
    }

    if let Some(description) = reference_sequence.description() {
        write_field(writer, tag::DESCRIPTION, description)?;
    }

    if let Some(md5_checksum) = reference_sequence.md5_checksum() {
        write_field(writer, tag::MD5_CHECKSUM, md5_checksum)?;
    }

    if let Some(species) = reference_sequence.species() {
        write_field(writer, tag::SPECIES, species)?;
    }

    if let Some(molecule_topology) = reference_sequence.molecule_topology() {
        write_field(writer, tag::MOLECULE_TOPOLOGY, molecule_topology)?;
    }

    if let Some(uri) = reference_sequence.uri() {
        write_field(writer, tag::URI, uri)?;
    }

    write_other_fields(writer, reference_sequence.other_fields())?;

    Ok(())
}
