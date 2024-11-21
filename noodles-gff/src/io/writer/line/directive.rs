use std::io::{self, Write};

use crate::{directive_buf::name, DirectiveBuf};

pub(super) fn write_directive<W>(writer: &mut W, directive: &DirectiveBuf) -> io::Result<()>
where
    W: Write,
{
    write_prefix(writer)?;

    match directive {
        DirectiveBuf::GffVersion(version) => {
            write_key(writer, name::GFF_VERSION)?;
            write_separator(writer)?;
            write!(writer, "{version}")?;
        }
        DirectiveBuf::SequenceRegion(sequence_region) => {
            write_key(writer, name::SEQUENCE_REGION)?;
            write_separator(writer)?;
            write!(writer, "{sequence_region}")?
        }
        DirectiveBuf::FeatureOntology(uri) => {
            write_key(writer, name::FEATURE_ONTOLOGY)?;
            write_separator(writer)?;
            write_value(writer, uri)?;
        }
        DirectiveBuf::AttributeOntology(uri) => {
            write_key(writer, name::ATTRIBUTE_ONTOLOGY)?;
            write_separator(writer)?;
            write_value(writer, uri)?;
        }
        DirectiveBuf::SourceOntology(uri) => {
            write_key(writer, name::SOURCE_ONTOLOGY)?;
            write_separator(writer)?;
            write_value(writer, uri)?;
        }
        DirectiveBuf::Species(uri) => {
            write_key(writer, name::SPECIES)?;
            write_separator(writer)?;
            write_value(writer, uri)?;
        }
        DirectiveBuf::GenomeBuild(genome_build) => {
            write_key(writer, name::GENOME_BUILD)?;
            write_separator(writer)?;
            write!(writer, "{genome_build}")?
        }
        DirectiveBuf::ForwardReferencesAreResolved => {
            write_key(writer, name::FORWARD_REFERENCES_ARE_RESOLVED)?
        }
        DirectiveBuf::StartOfFasta => write_key(writer, name::START_OF_FASTA)?,
        DirectiveBuf::Other(key, value) => {
            write_key(writer, key.as_ref())?;

            if let Some(v) = value {
                write_separator(writer)?;
                write_value(writer, v)?;
            }
        }
    }

    Ok(())
}

fn write_prefix<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const PREFIX: &[u8; 2] = b"##";
    writer.write_all(PREFIX)
}

fn write_key<W>(writer: &mut W, key: &str) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(key.as_bytes())
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: u8 = b' ';
    writer.write_all(&[SEPARATOR])
}

fn write_value<W>(writer: &mut W, value: &str) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(value.as_bytes())
}
