use std::io::{self, BufRead};

use super::DEFINITION_PREFIX;

pub(super) fn read_sequence<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    read_sequence_limit(reader, usize::MAX, buf)
}

pub(super) fn read_sequence_limit<R>(
    reader: &mut R,
    max_bases: usize,
    buf: &mut Vec<u8>,
) -> io::Result<usize>
where
    R: BufRead,
{
    use memchr::memchr;

    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    let mut len = 0;

    while buf.len() < max_bases {
        let src = reader.fill_buf()?;

        let is_eof = src.is_empty();
        let is_end_of_sequence = || src[0] == DEFINITION_PREFIX;

        if is_eof || is_end_of_sequence() {
            break;
        }

        let remaining_bases = max_bases - buf.len();

        let n = match memchr(LINE_FEED, src) {
            Some(i) => {
                let i = i.min(remaining_bases);
                let line = &src[..i];

                if line.ends_with(&[CARRIAGE_RETURN]) {
                    let end = line.len() - 1;
                    buf.extend(&line[..end]);
                } else {
                    buf.extend(line);
                }

                i + 1
            }
            None => {
                let i = remaining_bases.min(src.len());
                let line = &src[..i];
                buf.extend(line);
                i
            }
        };

        reader.consume(n);

        len += n;
    }

    Ok(len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_sequence() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, mut reader: &[u8], expected: &[u8]) -> io::Result<()> {
            buf.clear();
            read_sequence(&mut reader, buf)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, b"ACGT\n", b"ACGT")?;
        t(&mut buf, b"ACGT\n>sq1\n", b"ACGT")?;
        t(&mut buf, b"NNNN\nNNNN\nNN\n", b"NNNNNNNNNN")?;

        t(&mut buf, b"ACGT\r\n", b"ACGT")?;
        t(&mut buf, b"ACGT\r\n>sq1\r\n", b"ACGT")?;
        t(&mut buf, b"NNNN\r\nNNNN\r\nNN\r\n", b"NNNNNNNNNN")?;

        Ok(())
    }

    #[test]
    fn test_read_sequence_limit() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            mut reader: &[u8],
            max_bases: usize,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            read_sequence_limit(&mut reader, max_bases, buf)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, b"ACGT\n", 4, b"ACGT")?;
        t(&mut buf, b"ACGT\n>sq0\n", 4, b"ACGT")?;
        t(&mut buf, b"ACGT\nACGT\nAC\n", 10, b"ACGTACGTAC")?;

        t(&mut buf, b"ACGT\n", 2, b"AC")?;
        t(&mut buf, b"ACGT\n>sq0\n", 2, b"AC")?;
        t(&mut buf, b"ACGT\nACGT\nAC", 2, b"AC")?;

        t(&mut buf, b"ACGT\n", 5, b"ACGT")?;
        t(&mut buf, b"ACGT\n>sq0\n", 5, b"ACGT")?;
        t(&mut buf, b"ACGT\nACGT\nAC", 5, b"ACGTA")?;

        t(&mut buf, b"ACGT\n", 5, b"ACGT")?;
        t(&mut buf, b"ACGT\r\n>sq0\r\n", 5, b"ACGT")?;
        t(&mut buf, b"ACGT\r\nACGT\r\nAC\r\n", 5, b"ACGTA")?;

        Ok(())
    }
}
