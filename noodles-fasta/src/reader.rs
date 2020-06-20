use std::io::{self, BufRead, Seek, SeekFrom};

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

    pub fn read_definition(&mut self, buf: &mut String) -> io::Result<usize> {
        read_line(&mut self.inner, buf)
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

fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    let result = reader.read_line(buf);
    buf.pop();
    result
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn test_read_definition() -> io::Result<()> {
        let data = b">sq0\nACGT\n";
        let mut reader = Reader::new(&data[..]);

        let mut description_buf = String::new();
        reader.read_definition(&mut description_buf)?;

        assert_eq!(description_buf, ">sq0");

        Ok(())
    }

    #[test]
    fn test_read_sequence() -> io::Result<()> {
        let mut sequence_buf = Vec::new();

        let data = b"ACGT\n";
        let mut reader = Reader::new(&data[..]);
        reader.read_sequence(&mut sequence_buf)?;
        assert_eq!(sequence_buf, b"ACGT");

        let data = b"ACGT\n>sq1\n";
        let mut reader = Reader::new(&data[..]);
        sequence_buf.clear();
        reader.read_sequence(&mut sequence_buf)?;
        assert_eq!(sequence_buf, b"ACGT");

        let data = b"NNNN\nNNNN\nNN\n";
        let mut reader = Reader::new(&data[..]);
        sequence_buf.clear();
        reader.read_sequence(&mut sequence_buf)?;
        assert_eq!(sequence_buf, b"NNNNNNNNNN");

        Ok(())
    }

    #[test]
    fn test_read_sequence_after_seek() {
        let data = b">sq0\nACGT\n>sq1\nNNNN\n";
        let cursor = Cursor::new(&data[..]);
        let mut reader = Reader::new(cursor);

        reader.seek(SeekFrom::Start(14)).unwrap();

        let mut buf = Vec::new();
        reader.read_sequence(&mut buf).unwrap();

        assert_eq!(buf, b"NNNN");
    }
}
