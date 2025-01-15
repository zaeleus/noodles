use std::io::{self, BufRead, BufReader, Read};

use flate2::read::GzDecoder;

use crate::crai::Index;

/// A CRAM index reader.
pub struct Reader<R> {
    inner: BufReader<GzDecoder<R>>,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a CRAM index reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::crai;
    /// let data = [];
    /// let reader = crai::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner: BufReader::new(GzDecoder::new(inner)),
        }
    }

    /// Reads a CRAM index.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram::crai;
    /// let mut reader = File::open("sample.cram.crai").map(crai::io::Reader::new)?;
    /// let index = reader.read_index()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_index(&mut self) -> io::Result<Index> {
        let mut buf = String::new();
        let mut index = Vec::new();

        loop {
            buf.clear();

            match read_line(&mut self.inner, &mut buf) {
                Ok(0) => break,
                Ok(_) => {
                    let record = buf
                        .parse()
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                    index.push(record);
                }
                Err(e) => return Err(e),
            }
        }

        Ok(index)
    }
}

fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    match reader.read_line(buf) {
        Ok(0) => Ok(0),
        Ok(n) => {
            buf.pop();
            Ok(n)
        }
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use flate2::write::GzEncoder;
    use noodles_core::Position;

    use crate::crai::Record;

    use super::*;

    #[test]
    fn test_read_index() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"\
0\t10946\t6765\t17711\t233\t317811
0\t17711\t121393\t317811\t233\t317811
";

        let mut writer = GzEncoder::new(Vec::new(), Default::default());
        writer.write_all(data)?;
        let compressed_data = writer.finish()?;

        let mut reader = Reader::new(&compressed_data[..]);

        let actual = reader.read_index()?;

        let expected = vec![
            Record::new(Some(0), Position::new(10946), 6765, 17711, 233, 317811),
            Record::new(Some(0), Position::new(17711), 121393, 317811, 233, 317811),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
