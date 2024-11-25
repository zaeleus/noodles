mod value;

use std::io::{self, Write};

use self::value::write_value;
use crate::{directive_buf::key, DirectiveBuf};

pub(crate) fn write_directive<W>(writer: &mut W, directive: &DirectiveBuf) -> io::Result<()>
where
    W: Write,
{
    write_prefix(writer)?;

    match directive {
        DirectiveBuf::GffVersion(version) => {
            write_key(writer, key::GFF_VERSION)?;
            write_separator(writer)?;
            value::write_gff_version(writer, version)?;
        }
        DirectiveBuf::SequenceRegion(sequence_region) => {
            write_key(writer, key::SEQUENCE_REGION)?;
            write_separator(writer)?;
            value::write_sequence_region(writer, sequence_region)?;
        }
        DirectiveBuf::GenomeBuild(genome_build) => {
            write_key(writer, key::GENOME_BUILD)?;
            write_separator(writer)?;
            value::write_genome_build(writer, genome_build)?;
        }
        DirectiveBuf::Other(key, value) => {
            write_key(writer, key)?;

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
