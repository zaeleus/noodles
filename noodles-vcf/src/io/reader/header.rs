use std::io::{self, BufRead};

use crate::{header, Header};

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: BufRead,
{
    let mut parser = header::Parser::default();
    let mut buf = Vec::new();

    while read_header_line(reader, &mut buf)? != 0 {
        parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    parser
        .finish()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_header_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    const PREFIX: u8 = b'#';
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    let src = reader.fill_buf()?;

    if src.is_empty() || src[0] != PREFIX {
        return Ok(0);
    }

    dst.clear();

    match reader.read_until(LINE_FEED, dst)? {
        0 => Ok(0),
        n => {
            if dst.ends_with(&[LINE_FEED]) {
                dst.pop();

                if dst.ends_with(&[CARRIAGE_RETURN]) {
                    dst.pop();
                }
            }

            Ok(n)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn collect_lines<R>(reader: &mut R) -> io::Result<Vec<Vec<u8>>>
    where
        R: BufRead,
    {
        let mut buf = Vec::new();
        let mut lines = Vec::new();

        while read_header_line(reader, &mut buf)? != 0 {
            lines.push(buf.clone());
        }

        Ok(lines)
    }

    #[test]
    fn test_read_raw_header() -> io::Result<()> {
        let src = b"##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
sq0\t1\t.\tA\t.\t.\tPASS\t.
";

        let mut reader = &src[..];

        let actual = collect_lines(&mut reader)?;
        let expected = [
            b"##fileformat=VCFv4.3".to_vec(),
            b"##fileDate=20200501".to_vec(),
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_vec(),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_raw_header_with_no_records() -> io::Result<()> {
        let src = b"##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        let mut reader = &src[..];

        let actual = collect_lines(&mut reader)?;
        let expected = [
            b"##fileformat=VCFv4.3".to_vec(),
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_vec(),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_raw_header_with_multiple_buffer_fills() -> io::Result<()> {
        use std::io::BufReader;

        let src = b"##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        let mut reader = BufReader::with_capacity(16, &src[..]);
        let actual = collect_lines(&mut reader)?;
        let expected = [
            b"##fileformat=VCFv4.3".to_vec(),
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_vec(),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_raw_header_with_no_header() -> io::Result<()> {
        let src = [];
        let mut reader = &src[..];
        let actual = collect_lines(&mut reader)?;
        assert!(actual.is_empty());

        let src = b"sq0\t1\t.\tA\t.\t.\tPASS\t.\n";
        let mut reader = &src[..];
        let actual = collect_lines(&mut reader)?;
        assert!(actual.is_empty());

        Ok(())
    }

    #[test]
    fn test_read_raw_header_with_missing_end_of_line() -> io::Result<()> {
        let src = b"##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

        let mut reader = &src[..];

        let actual = collect_lines(&mut reader)?;
        let expected = [
            b"##fileformat=VCFv4.3".to_vec(),
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_vec(),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
