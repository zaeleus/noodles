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

#[cfg(test)]
mod tests {
    use crate::directive_buf::value::{GenomeBuild, SequenceRegion};

    use super::*;

    #[test]
    fn test_write_directive() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, directive: &DirectiveBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_directive(buf, directive)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(
            &mut buf,
            &DirectiveBuf::new(
                key::GFF_VERSION,
                Some(Value::GffVersion(Default::default())),
            ),
            b"##gff-version 3",
        )?;

        t(
            &mut buf,
            &DirectiveBuf::new(
                key::SEQUENCE_REGION,
                Some(Value::SequenceRegion(SequenceRegion::new(
                    String::from("sq0"),
                    8,
                    13,
                ))),
            ),
            b"##sequence-region sq0 8 13",
        )?;

        t(
            &mut buf,
            &DirectiveBuf::new(
                key::GENOME_BUILD,
                Some(Value::GenomeBuild(GenomeBuild::new(
                    String::from("NDLS"),
                    String::from("r1"),
                ))),
            ),
            b"##genome-build NDLS r1",
        )?;

        t(
            &mut buf,
            &DirectiveBuf::new("noodles", Some(Value::from("gff"))),
            b"##noodles gff",
        )?;

        t(&mut buf, &DirectiveBuf::new("noodles", None), b"##noodles")?;

        buf.clear();
        assert!(matches!(
            write_directive(
                &mut buf,
                &DirectiveBuf::new(
                    key::GENOME_BUILD,
                    Some(Value::GffVersion(Default::default())),
                ),
            ),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
