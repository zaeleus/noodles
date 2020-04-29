use std::io::{self, BufRead};

use super::Record;

pub struct Reader<R> {
    inner: R,
    line_buf: String,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            line_buf: String::new(),
        }
    }

    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        self.line_buf.clear();

        match self.inner.read_line(&mut self.line_buf) {
            Ok(0) => Ok(0),
            Ok(n) => {
                self.line_buf.pop();

                *record = self
                    .line_buf
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
sq0\t10946\t4\t80\t81
sq1\t17711\t10954\t80\t81
";

        let mut reader = Reader::new(&data[..]);
        let mut record = Record::default();

        let bytes_read = reader.read_record(&mut record)?;
        assert_eq!(bytes_read, 18);

        let bytes_read = reader.read_record(&mut record)?;
        assert_eq!(bytes_read, 22);

        let bytes_read = reader.read_record(&mut record)?;
        assert_eq!(bytes_read, 0);

        Ok(())
    }
}
