use std::io::{self, BufRead};

use crate::{header, Header};

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: BufRead,
{
    let mut parser = header::Parser::default();
    let mut buf = String::new();

    while read_header_line(reader, &mut buf)? != 0 {
        parser
            .parse_partial(buf.as_bytes())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    Ok(parser.finish())
}

fn read_header_line<R>(reader: &mut R, dst: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    const PREFIX: u8 = b'@';
    const LINE_FEED: char = '\n';
    const CARRIAGE_RETURN: char = '\r';

    let src = reader.fill_buf()?;

    if src.is_empty() || src[0] != PREFIX {
        return Ok(0);
    }

    dst.clear();

    match reader.read_line(dst)? {
        0 => Ok(0),
        n => {
            if dst.ends_with(LINE_FEED) {
                dst.pop();

                if dst.ends_with(CARRIAGE_RETURN) {
                    dst.pop();
                }
            }

            Ok(n)
        }
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use super::*;
    use crate::header::record::value::{
        map::{self, header::Version, ReferenceSequence},
        Map,
    };

    #[test]
    fn test_read_header_with_no_header() -> io::Result<()> {
        let data = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n";
        let mut reader = &data[..];
        assert!(read_header(&mut reader)?.is_empty());
        Ok(())
    }

    #[test]
    fn test_read_header_with_no_records() -> io::Result<()> {
        let data = "@HD\tVN:1.6\n";
        let mut reader = data.as_bytes();

        let actual = read_header(&mut reader)?;

        let expected = crate::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_with_multiple_buffer_fills() -> Result<(), Box<dyn std::error::Error>> {
        use std::io::BufReader;

        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(ln) => ln,
            None => unreachable!(),
        };

        let data = "@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:8\n";
        let mut reader = BufReader::with_capacity(16, data.as_bytes());

        let actual = read_header(&mut reader)?;

        let expected = crate::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence("sq0".parse()?, Map::<ReferenceSequence>::new(SQ0_LN))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
