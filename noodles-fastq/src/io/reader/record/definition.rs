use std::io::{self, BufRead};

use super::{CARRIAGE_RETURN, LINE_FEED, read_line, read_u8};
use crate::record::Definition;

pub(crate) fn read_definition<R>(reader: &mut R, definition: &mut Definition) -> io::Result<usize>
where
    R: BufRead,
{
    use memchr::memchr3;

    const NAME_PREFIX: u8 = b'@';

    const HORIZONTAL_TAB: u8 = b'\t';
    const SPACE: u8 = b' ';

    match read_u8(reader) {
        Ok(prefix) => {
            if prefix != NAME_PREFIX {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid name prefix",
                ));
            }
        }
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
        Err(e) => return Err(e),
    }

    let mut is_eol = false;
    let mut len = 1;

    loop {
        let src = reader.fill_buf()?;

        if src.is_empty() {
            break;
        }

        let (matched_needle, n) = match memchr3(SPACE, HORIZONTAL_TAB, LINE_FEED, src) {
            Some(i) => {
                let name_src = match src[i] {
                    SPACE | HORIZONTAL_TAB => &src[..i],
                    LINE_FEED => {
                        is_eol = true;

                        let line = &src[..i];

                        if line.ends_with(&[CARRIAGE_RETURN]) {
                            // SAFETY: `line.len()` is > 0.
                            let end = line.len() - 1;
                            &line[..end]
                        } else {
                            line
                        }
                    }
                    _ => unreachable!(),
                };

                definition.name_mut().extend(name_src);

                (true, i + 1)
            }
            None => {
                definition.name_mut().extend(src);
                (false, src.len())
            }
        };

        len += n;
        reader.consume(n);

        if matched_needle {
            break;
        }
    }

    if !is_eol {
        len += read_line(reader, definition.description_mut())?;
    }

    Ok(len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_definition() -> io::Result<()> {
        let mut definition = Definition::default();

        let data = b"@r0\n";
        let mut reader = &data[..];
        definition.clear();
        read_definition(&mut reader, &mut definition)?;
        assert_eq!(definition.name(), &b"r0"[..]);
        assert!(definition.description().is_empty());

        let data = b"@r0 LN:4\n";
        let mut reader = &data[..];
        definition.clear();
        read_definition(&mut reader, &mut definition)?;
        assert_eq!(definition.name(), &b"r0"[..]);
        assert_eq!(definition.description(), &b"LN:4"[..]);

        let data = b"@r0\tLN:4\n";
        let mut reader = &data[..];
        definition.clear();
        read_definition(&mut reader, &mut definition)?;
        assert_eq!(definition.name(), &b"r0"[..]);
        assert_eq!(definition.description(), &b"LN:4"[..]);

        // https://github.com/zaeleus/noodles/issues/166
        let data = b"@\nA\r";
        let mut reader = &data[..];
        definition.clear();
        read_definition(&mut reader, &mut definition)?;
        assert!(definition.name().is_empty());
        assert!(definition.description().is_empty());

        let data = b"r0\n";
        let mut reader = &data[..];
        definition.clear();
        assert!(matches!(
            read_definition(&mut reader, &mut definition),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
