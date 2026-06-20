use std::io::{self, BufRead};

use super::read_line;
use crate::record::Definition;

pub(super) fn read_definition<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    definition: &mut Definition,
) -> io::Result<usize>
where
    R: BufRead,
{
    buf.clear();

    match read_line(reader, buf)? {
        0 => Ok(0),
        n => {
            let (name, description) = parse_definition(buf)?;

            let description = if description.is_empty() {
                None
            } else {
                Some(description.into())
            };

            *definition = Definition::new(name, description);

            Ok(n)
        }
    }
}

pub(crate) fn parse_definition(mut src: &[u8]) -> io::Result<(&[u8], &[u8])> {
    if src.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "empty input"));
    }

    consume_prefix(&mut src)?;
    let name = parse_name(&mut src)?;
    let description = src.trim_ascii();

    Ok((name, description))
}

fn consume_prefix(src: &mut &[u8]) -> io::Result<()> {
    const PREFIX: u8 = b'>';

    let b = src
        .split_off_first()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing prefix"))?;

    if *b == PREFIX {
        Ok(())
    } else {
        Err(io::Error::new(io::ErrorKind::InvalidData, "invalid prefix"))
    }
}

fn parse_name<'a>(src: &mut &'a [u8]) -> io::Result<&'a [u8]> {
    let i = src
        .iter()
        .position(|&b| b.is_ascii_whitespace())
        .unwrap_or(src.len());

    // SAFETY: `i <= src.len()`.
    let name = src.split_off(..i).unwrap();

    if name.is_empty() {
        Err(io::Error::new(io::ErrorKind::InvalidData, "missing name"))
    } else {
        Ok(name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_definition() -> io::Result<()> {
        assert_eq!(parse_definition(&b">sq0"[..])?, (&b"sq0"[..], &b""[..]));
        assert_eq!(
            parse_definition(&b">sq0  LN:13"[..])?,
            (&b"sq0"[..], &b"LN:13"[..])
        );

        assert!(matches!(
            parse_definition(b""),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
        assert!(matches!(
            parse_definition(b"sq0"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
        assert!(matches!(
            parse_definition(b">"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
