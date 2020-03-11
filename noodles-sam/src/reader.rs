use std::io::{self, BufRead};

#[derive(Debug)]
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    pub fn read_header(&mut self) -> io::Result<Vec<u8>> {
        let mut header_buf = Vec::new();
        let mut eol = false;

        loop {
            let buf = self.inner.fill_buf()?;

            if eol && buf.len() > 0 && buf[0] != b'@' {
                break;
            }

            let (read_eol, len) = match buf.iter().position(|&b| b == b'\n') {
                Some(i) => {
                    header_buf.extend(&buf[..=i]);
                    (true, i + 1)
                }
                None => {
                    header_buf.extend(buf);
                    (false, buf.len())
                }
            };

            eol = read_eol;
            self.inner.consume(len);
        }

        Ok(header_buf)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static SAM: &[u8] = b"\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:ref0\tLN:13
r001
";

    #[test]
    fn test_read_header() -> io::Result<()> {
        let mut reader = Reader::new(&SAM[..]);

        let actual = reader.read_header()?;
        let expected = b"\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:ref0\tLN:13
";
        assert_eq!(actual, &expected[..]);

        Ok(())
    }
}
