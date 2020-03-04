use std::{
    io::{self, BufRead},
    str,
};

use crate::Record;

pub struct Reader<R> {
    inner: R,
    buf: Vec<u8>,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            buf: Vec::new(),
        }
    }

    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        record.clear();
        *record.name_mut() = self.read_name()?;
        self.read_sequence(record.sequence_mut())
    }

    pub fn read_name(&mut self) -> io::Result<String> {
        if self.buf.is_empty() {
            match read_line(&mut self.inner, &mut self.buf) {
                Ok(0) => return Ok(String::new()),
                Ok(_) => {}
                Err(e) => return Err(e),
            }
        }

        let description = str::from_utf8(&self.buf[1..])
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        description
            .split_whitespace()
            .next()
            .map(|name| name.into())
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing sequence name"))
    }

    pub fn read_sequence(&mut self, buf: &mut Vec<u8>) -> io::Result<usize> {
        let mut bytes_read = 0;

        loop {
            match read_line(&mut self.inner, &mut self.buf) {
                Ok(0) => return Ok(bytes_read),
                Ok(n) => bytes_read += n,
                Err(e) => return Err(e),
            }

            if self.buf[0] == b'>' {
                bytes_read -= buf.len();
                break;
            }

            buf.extend(&self.buf);
        }

        Ok(bytes_read)
    }
}

fn read_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    buf.clear();
    let result = reader.read_until(b'\n', buf);
    buf.pop();
    result
}

#[cfg(test)]
mod tests {
    use crate::Record;

    use super::*;

    #[test]
    fn test_read_sequence() {
        let data = b"\
>chr1
ATCG
>chr2
NNNN
NNNNNN
";

        let mut reader = Reader::new(&data[..]);
        let mut record = Record::default();

        reader.read_record(&mut record).unwrap();
        assert_eq!(record.name(), "chr1");
        assert_eq!(record.sequence(), b"ATCG");

        reader.read_record(&mut record).unwrap();
        assert_eq!(record.name(), "chr2");
        assert_eq!(record.sequence(), b"NNNNNNNNNN");
    }

    #[test]
    fn test_read_sequence_with_empty_description() {
        let data = b">\nATCG";
        let mut reader = Reader::new(&data[..]);
        let mut record = Record::default();
        assert!(reader.read_record(&mut record).is_err());
    }
}
