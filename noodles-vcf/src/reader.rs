use std::io::{self, BufRead};

const NEWLINE: u8 = b'\n';
const META_PREFIX: &[u8] = b"##";
const HEADER_PREFIX: char = '#';

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

    pub fn read_meta(&mut self) -> io::Result<String> {
        let mut meta_buf = Vec::new();
        let mut eol = false;

        loop {
            let buf = self.inner.fill_buf()?;

            if eol && buf.len() > 2 && &buf[..2] != META_PREFIX {
                break;
            }

            let (read_eol, len) = match buf.iter().position(|&b| b == NEWLINE) {
                Some(i) => {
                    meta_buf.extend(&buf[..=i]);
                    (true, i + 1)
                }
                None => {
                    meta_buf.extend(buf);
                    (false, buf.len())
                }
            };

            eol = read_eol;
            self.inner.consume(len);
        }

        String::from_utf8(meta_buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    pub fn read_header(&mut self) -> io::Result<String> {
        let mut buf = String::new();
        self.inner.read_line(&mut buf)?;

        if !buf.starts_with(HEADER_PREFIX) {
            Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "header missing prefix",
            ))
        } else {
            buf.pop();
            Ok(buf)
        }
    }

    pub fn read_record(&mut self, buf: &mut String) -> io::Result<usize> {
        let result = self.inner.read_line(buf);
        buf.pop();
        result
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
    fn test_read_meta() -> io::Result<()> {
        let mut reader = Reader::new(DATA);

        let actual = reader.read_meta()?;
        let expected = "\
##fileformat=VCFv4.3
##fileDate=20200501
";

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header() -> io::Result<()> {
        let mut reader = Reader::new(DATA);
        reader.read_meta()?;

        let actual = reader.read_header()?;
        let expected = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_with_invalid_prefix() {
        let data = b"CHROM";
        let mut reader = Reader::new(&data[..]);
        assert!(reader.read_header().is_err());
    }

    #[test]
    fn test_read_record() -> io::Result<()> {
        let mut reader = Reader::new(DATA);
        reader.read_meta()?;
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
