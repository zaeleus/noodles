use std::io::{self, BufRead};

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<String>
where
    R: BufRead,
{
    use memchr::memchr;

    const PREFIX: u8 = b'@';
    const LINE_FEED: u8 = b'\n';

    let mut header_buf = Vec::new();
    let mut is_eol = false;

    for i in 0.. {
        let buf = reader.fill_buf()?;

        if (i == 0 || is_eol) && buf.first().map(|&b| b != PREFIX).unwrap_or(true) {
            break;
        }

        let (read_eol, len) = if let Some(i) = memchr(LINE_FEED, buf) {
            header_buf.extend(&buf[..=i]);
            (true, i + 1)
        } else {
            header_buf.extend(buf);
            (false, buf.len())
        };

        is_eol = read_eol;

        reader.consume(len);
    }

    String::from_utf8(header_buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let header = read_header(&mut reader)?;
        assert_eq!(header, data);
        Ok(())
    }

    #[test]
    fn test_read_header_with_multiple_buffer_fills() -> io::Result<()> {
        use std::io::BufReader;

        let data = "@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:8\n";
        let mut reader = BufReader::with_capacity(16, data.as_bytes());
        let header = read_header(&mut reader)?;

        assert_eq!(header, data);

        Ok(())
    }
}
