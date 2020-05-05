use std::io::{self, BufRead};

use super::Record;

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

    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        let mut buf = String::new();

        match self.inner.read_line(&mut buf) {
            Ok(0) => Ok(0),
            Ok(n) => {
                buf.pop();

                *record = buf
                    .parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                Ok(n)
            }
            Err(e) => Err(e),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_record() -> io::Result<()> {
        let data = b"\
0\t10946\t6765\t17711\t233\t317811
0\t17711\t121393\t317811\t233\t317811
";

        let mut reader = Reader::new(&data[..]);
        let mut record = Record::default();

        let bytes_read = reader.read_record(&mut record)?;
        assert_eq!(bytes_read, 30);

        let bytes_read = reader.read_record(&mut record)?;
        assert_eq!(bytes_read, 33);

        let bytes_read = reader.read_record(&mut record)?;
        assert_eq!(bytes_read, 0);

        Ok(())
    }
}
