mod definition;

pub(crate) use self::definition::read_definition;

use std::io::{self, BufRead, Read};

use crate::Record;

const LINE_FEED: u8 = b'\n';
const CARRIAGE_RETURN: u8 = b'\r';

pub(super) fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: BufRead,
{
    record.clear();

    let mut len = match read_definition(reader, record.definition_mut()) {
        Ok(0) => return Ok(0),
        Ok(n) => n,
        Err(e) => return Err(e),
    };

    len += read_line(reader, record.sequence_mut())?;
    len += consume_plus_line(reader)?;
    len += read_line(reader, record.quality_scores_mut())?;

    Ok(len)
}

fn read_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    match reader.read_until(LINE_FEED, buf) {
        Ok(0) => Ok(0),
        Ok(n) => {
            if buf.ends_with(&[LINE_FEED]) {
                buf.pop();

                if buf.ends_with(&[CARRIAGE_RETURN]) {
                    buf.pop();
                }
            }

            Ok(n)
        }
        Err(e) => Err(e),
    }
}

fn consume_line<R>(reader: &mut R) -> io::Result<usize>
where
    R: BufRead,
{
    use memchr::memchr;

    let mut is_eol = false;
    let mut len = 0;

    loop {
        let src = reader.fill_buf()?;

        if src.is_empty() || is_eol {
            break;
        }

        let n = match memchr(LINE_FEED, src) {
            Some(i) => {
                is_eol = true;
                i + 1
            }
            None => src.len(),
        };

        reader.consume(n);

        len += n;
    }

    Ok(len)
}

fn read_u8<R>(reader: &mut R) -> io::Result<u8>
where
    R: Read,
{
    let mut buf = [0; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn consume_plus_line<R>(reader: &mut R) -> io::Result<usize>
where
    R: BufRead,
{
    const PREFIX: u8 = b'+';

    match read_u8(reader)? {
        PREFIX => consume_line(reader).map(|n| n + 1),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid description prefix",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_line() -> io::Result<()> {
        let mut buf = Vec::new();

        let data = b"noodles\n";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, b"noodles");

        let data = b"noodles\r\n";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, b"noodles");

        let data = b"noodles";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, b"noodles");

        Ok(())
    }

    #[test]
    fn test_consume_plus_line() -> io::Result<()> {
        let data = b"+r0\n";
        let mut reader = &data[..];
        consume_plus_line(&mut reader)?;

        let data = b"r0\n";
        let mut reader = &data[..];
        assert!(matches!(
            consume_plus_line(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
