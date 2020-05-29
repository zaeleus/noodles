mod records;

pub use self::records::Records;

use std::io::{self, BufRead};

const NEWLINE: u8 = b'\n';
const HEADER_PREFIX: u8 = b'#';

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

    pub fn read_header(&mut self) -> io::Result<String> {
        let mut header_buf = Vec::new();
        let mut eol = false;

        loop {
            let buf = match self.inner.fill_buf() {
                Ok(buf) => buf,
                Err(ref e) if e.kind() == io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e),
            };

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

        String::from_utf8(header_buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
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
##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
sq0\t8
sq0\t13
";

    #[test]
    fn test_read_header() -> io::Result<()> {
        let mut reader = Reader::new(DATA);

        let actual = reader.read_header()?;
        let expected = "\
##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_record() -> io::Result<()> {
        let mut reader = Reader::new(DATA);
        reader.read_header()?;

        let mut buf = String::new();
        reader.read_record(&mut buf)?;
        assert_eq!(buf, "sq0\t8");

        buf.clear();
        reader.read_record(&mut buf)?;
        assert_eq!(buf, "sq0\t13");

        buf.clear();
        let len = reader.read_record(&mut buf)?;
        assert_eq!(len, 0);

        Ok(())
    }
}
