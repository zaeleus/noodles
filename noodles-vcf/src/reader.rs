use std::io::{self, BufRead};

const NEWLINE: u8 = b'\n';
const META_PREFIX: &[u8] = b"##";

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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_meta() -> io::Result<()> {
        let data = b"\
##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        let mut reader = Reader::new(&data[..]);

        let actual = reader.read_meta()?;
        let expected = "\
##fileformat=VCFv4.3
##fileDate=20200501
";

        assert_eq!(actual, expected);

        Ok(())
    }
}
