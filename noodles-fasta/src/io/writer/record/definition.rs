use std::io::{self, Write};

use crate::record::Definition;

const SEPARATOR: u8 = b' ';
const PREFIX: u8 = b'>';

pub(super) fn write_definition<W>(writer: &mut W, definition: &Definition) -> io::Result<()>
where
    W: Write,
{
    write_prefix(writer)?;
    write_name(writer, definition.name())?;

    if let Some(description) = definition.description() {
        write_separator(writer)?;
        write_description(writer, description)?;
    }

    Ok(())
}

fn write_prefix<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[PREFIX])
}

fn write_name<W>(writer: &mut W, name: &[u8]) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(name)
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[SEPARATOR])
}

fn write_description<W>(writer: &mut W, description: &[u8]) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(description)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_definition() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, definition: &Definition, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_definition(buf, definition)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Definition::new("sq0", None), b">sq0")?;
        t(
            &mut buf,
            &Definition::new("sq0", Some(Vec::from("LN:8"))),
            b">sq0 LN:8",
        )?;

        Ok(())
    }
}
