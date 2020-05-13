mod records;

pub use self::records::Records;

use std::io::{self, BufRead};

const HEADER_PREFIX: u8 = b'@';
const NEWLINE: u8 = b'\n';

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

            if eol && !buf.is_empty() && buf[0] != HEADER_PREFIX {
                break;
            }

            let (read_eol, len) = match buf.iter().position(|&b| b == NEWLINE) {
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

    pub fn read_record(&mut self, buf: &mut String) -> io::Result<usize> {
        let result = self.inner.read_line(buf);
        buf.pop();
        result
    }

    pub fn records(&mut self) -> Records<R> {
        Records::new(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static DATA: &[u8] = b"\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:ref0\tLN:13
r001
r002
";

    #[test]
    fn test_read_header() -> io::Result<()> {
        let mut reader = Reader::new(&DATA[..]);

        let actual = reader.read_header()?;
        let expected = b"\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:ref0\tLN:13
";
        assert_eq!(actual, &expected[..]);

        Ok(())
    }

    #[test]
    fn test_read_record() -> io::Result<()> {
        let mut reader = Reader::new(&DATA[..]);
        reader.read_header()?;

        let mut buf = String::new();
        reader.read_record(&mut buf)?;
        assert_eq!(buf, "r001");

        buf.clear();
        reader.read_record(&mut buf)?;
        assert_eq!(buf, "r002");

        buf.clear();
        let len = reader.read_record(&mut buf)?;
        assert_eq!(len, 0);

        Ok(())
    }
}
