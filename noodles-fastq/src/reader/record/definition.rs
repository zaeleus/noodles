use std::io::{self, BufRead};

use super::{read_line, read_u8, CARRIAGE_RETURN, LINE_FEED};
use crate::Record;

pub(crate) fn read_definition<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: BufRead,
{
    use memchr::memchr2;

    const DELIMITER: u8 = b' ';
    const NAME_PREFIX: u8 = b'@';

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

        let (matched_needle, n) = match memchr2(DELIMITER, LINE_FEED, src) {
            Some(i) => {
                let name_src = match src[i] {
                    DELIMITER => &src[..i],
                    LINE_FEED => {
                        is_eol = true;

                        if src.ends_with(&[CARRIAGE_RETURN]) {
                            &src[..i - 1]
                        } else {
                            &src[..i]
                        }
                    }
                    _ => unreachable!(),
                };

                record.name_mut().extend(name_src);

                (true, i + 1)
            }
            None => {
                record.name_mut().extend(src);
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
        len += read_line(reader, record.description_mut())?;
    }

    Ok(len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_definition() -> io::Result<()> {
        let mut record = Record::default();

        let data = b"@r0\n";
        let mut reader = &data[..];
        record.clear();
        read_definition(&mut reader, &mut record)?;
        assert_eq!(record.name(), b"r0");
        assert!(record.description().is_empty());

        let data = b"@r0 LN:4\n";
        let mut reader = &data[..];
        record.clear();
        read_definition(&mut reader, &mut record)?;
        assert_eq!(record.name(), b"r0");
        assert_eq!(record.description(), b"LN:4");

        let data = b"r0\n";
        let mut reader = &data[..];
        record.clear();
        assert!(matches!(
            read_definition(&mut reader, &mut record),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
