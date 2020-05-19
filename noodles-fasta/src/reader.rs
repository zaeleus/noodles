use std::io::{self, BufRead, Seek, SeekFrom};

use crate::Record;

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
        record.clear();

        let mut bytes_read = match self.read_description(record.name_mut()) {
            Ok(0) => return Ok(0),
            Ok(len) => len,
            Err(e) => return Err(e),
        };

        bytes_read += self.read_sequence(record.sequence_mut())?;

        Ok(bytes_read)
    }

    pub fn read_description(&mut self, buf: &mut String) -> io::Result<usize> {
        match self.inner.read_line(buf) {
            Ok(0) => Ok(0),
            Ok(len) => {
                if !buf.starts_with('>') {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "line does not start with a '>'",
                    ));
                }

                buf.pop();

                Ok(len)
            }
            Err(e) => Err(e),
        }
    }

    pub fn read_sequence(&mut self, buf: &mut Vec<u8>) -> io::Result<usize> {
        let mut bytes_read = 0;

        loop {
            let reader_buf = self.inner.fill_buf()?;

            if reader_buf.is_empty() || reader_buf[0] == b'>' {
                break;
            }

            let len = match reader_buf.iter().position(|&b| b == b'\n') {
                Some(i) => {
                    buf.extend(&reader_buf[..i]);
                    i + 1
                }
                None => {
                    buf.extend(reader_buf);
                    reader_buf.len()
                }
            };

            self.inner.consume(len);

            bytes_read += len;
        }

        Ok(bytes_read)
    }
}

impl<R> Seek for Reader<R>
where
    R: BufRead + Seek,
{
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.inner.seek(pos)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use crate::Record;

    use super::*;

    static DATA: &[u8] = b"\
>chr1
ATCG
>chr2
NNNN
NNNN
NN
";

    #[test]
    fn test_read_record() {
        let mut reader = Reader::new(&DATA[..]);
        let mut record = Record::default();

        reader.read_record(&mut record).unwrap();
        assert_eq!(record.name(), ">chr1");
        assert_eq!(record.sequence(), b"ATCG");

        reader.read_record(&mut record).unwrap();
        assert_eq!(record.name(), ">chr2");
        assert_eq!(record.sequence(), b"NNNNNNNNNN");
    }

    #[test]
    fn test_read_sequence_after_seek() {
        let cursor = Cursor::new(&DATA[..]);
        let mut reader = Reader::new(cursor);

        reader.seek(SeekFrom::Start(17)).unwrap();

        let mut buf = Vec::new();
        reader.read_sequence(&mut buf).unwrap();

        assert_eq!(buf, b"NNNNNNNNNN");
    }
}
