mod value;

use std::io::{self, Write};

use self::value::write_value;
use crate::{
    directive_buf::{key, Value},
    DirectiveBuf,
};

pub(crate) fn write_directive<W>(writer: &mut W, directive: &DirectiveBuf) -> io::Result<()>
where
    W: Write,
{
    write_prefix(writer)?;

    match (directive.key(), directive.value()) {
        (key::GFF_VERSION, Some(Value::GffVersion(version))) => {
            write_key(writer, key::GFF_VERSION)?;
            write_separator(writer)?;
            value::write_gff_version(writer, version)?;
        }
        (key::SEQUENCE_REGION, Some(Value::SequenceRegion(sequence_region))) => {
            write_key(writer, key::SEQUENCE_REGION)?;
            write_separator(writer)?;
            value::write_sequence_region(writer, sequence_region)?;
        }
        (key::GENOME_BUILD, Some(Value::GenomeBuild(genome_build))) => {
            write_key(writer, key::GENOME_BUILD)?;
            write_separator(writer)?;
            value::write_genome_build(writer, genome_build)?;
        }
        (key, Some(Value::String(value))) => {
            write_key(writer, key)?;
            write_separator(writer)?;
            write_value(writer, value)?;
        }
        (key, None) => {
            write_key(writer, key)?;
        }
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid directive",
            ))
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
