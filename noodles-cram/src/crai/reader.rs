use std::io::{self, BufRead};

use super::Index;

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
    use std::convert::TryFrom;

    use noodles_bam as bam;

    use crate::crai::Record;

    use super::*;

    #[test]
    fn test_read_index() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"\
0\t10946\t6765\t17711\t233\t317811
0\t17711\t121393\t317811\t233\t317811
";

        let mut reader = Reader::new(&data[..]);

        let actual = reader.read_index()?;

        let expected = vec![
            Record::new(
                bam::record::ReferenceSequenceId::try_from(0).map(Some)?,
                10946,
                6765,
                17711,
                233,
                317811,
            ),
            Record::new(
                bam::record::ReferenceSequenceId::try_from(0).map(Some)?,
                17711,
                121393,
                317811,
                233,
                317811,
            ),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
