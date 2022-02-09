//! FASTX record

/* std use */

/* crate use */
use bstr::io::BufReadExt;
use bstr::ByteSlice;

/* project use */
use crate::error::Error;
use crate::record::Record;
use crate::records::Records;

/// Struct to perform reade of fasta or fastq file
pub struct Reader<R>
where
    R: std::io::BufRead,
{
    inner: bstr::io::ByteLines<R>,
    previous_line: Option<Vec<u8>>,
}

impl<R> Reader<R>
where
    R: std::io::BufRead,
{
    /// Create a new Reader from a [std::io::BufRead]
    pub fn new(inner: R) -> std::io::Result<Self> {
        let mut iter = inner.byte_lines();
        let previous_line = iter
            .next()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::UnexpectedEof, ""))??;

        Ok(Self {
            inner: iter,
            previous_line: Some(previous_line),
        })
    }

    /// Get next record from inner
    pub fn next_record(&mut self) -> Option<std::io::Result<Record>> {
        match self.previous_line.clone() {
            Some(line) if line[0] == b'>' => self.read_fasta(&line),
            Some(line) if line[0] == b'@' => self.read_fastq(&line),
            Some(_) => Some(Err(Error::MissingPrefix.into())),
            None => None,
        }
    }

    /// Generate an iterator of Record
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self)
    }

    fn read_fasta(&mut self, previous_line: &[u8]) -> Option<std::io::Result<Record>> {
        let mut record = Record::default();

        if previous_line.len() <= 1 {
            return Some(Err(Error::MissingName.into()));
        }
        match previous_line.find_byte(b' ') {
            Some(i) => {
                record.name_mut().extend(&previous_line[1..i]);
                *record.description_mut() = Some(previous_line[i + 1..].to_vec());
            }
            None => record.name_mut().extend(&previous_line[1..]),
        }

        loop {
            let line = match self.inner.next() {
                Some(Ok(line)) => line,
                Some(Err(e)) => return Some(Err(e)),
                None => {
                    self.previous_line = None;
                    break;
                }
            };

            if line.starts_with(b"@") || line.starts_with(b">") {
                self.previous_line = Some(line);
                break;
            }

            record.sequence_mut().extend(&line);
        }

        if record.sequence().is_empty() {
            return Some(Err(Error::MissingSequence.into()));
        }

        Some(Ok(record))
    }

    fn read_fastq(&mut self, previous_line: &[u8]) -> Option<std::io::Result<Record>> {
        let mut record = Record::default();

        if previous_line.len() <= 1 {
            return Some(Err(Error::MissingName.into()));
        }

        // Read header
        match previous_line.find_byte(b' ') {
            Some(i) => {
                record.name_mut().extend(&previous_line[1..i]);
                *record.description_mut() = Some(previous_line[i + 1..].to_vec());
            }
            None => record.name_mut().extend(&previous_line[1..]),
        }

        // Read sequence
        let seq_line = match self.inner.next() {
            Some(Ok(a)) => a,
            Some(Err(e)) => return Some(Err(e)),
            None => return Some(Err(Error::PartialRecord.into())),
        };
        if seq_line.is_empty() {
            return Some(Err(Error::MissingSequence.into()));
        }
        *record.sequence_mut() = seq_line;

        // Second description
        let second_desc = match self.inner.next() {
            Some(Ok(a)) => a,
            Some(Err(e)) => return Some(Err(e)),
            None => return Some(Err(Error::PartialRecord.into())),
        };
        if second_desc[0] != b'+' {
            return Some(Err(Error::MissingSecondDescription.into()));
        }
        if second_desc.len() <= 1 {
            *record.second_description_mut() = None;
        } else {
            *record.second_description_mut() = Some(second_desc[1..].to_vec());
        }

        // quality
        let qual = match self.inner.next() {
            Some(Ok(a)) => a,
            Some(Err(e)) => return Some(Err(e)),
            None => return Some(Err(Error::PartialRecord.into())),
        };
        if qual.is_empty() {
            return Some(Err(Error::MissingQuality.into()));
        }
        *record.quality_mut() = Some(qual[..].to_vec());

        match self.inner.next() {
            Some(Ok(line)) => self.previous_line = Some(line),
            Some(Err(e)) => return Some(Err(e)),
            None => self.previous_line = None,
        }

        Some(Ok(record))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_record() -> std::io::Result<()> {
        let data = b"\
@noodles:1/1
AGCT
+
abcd
@noodles: 2/1
TCGA
+noodles:2/1
dcba
>fasta sequence in middle of my fastq
GATC
ATGA
>1
ACGT
";

        let mut reader = Reader::new(&data[..])?;

        let record = reader.next_record().unwrap()?;
        assert_eq!(b"noodles:1/1", record.name());
        assert_eq!(None, record.description());
        assert_eq!(b"AGCT", record.sequence());
        assert_eq!(None, record.second_description());
        assert_eq!(Some(b"abcd".as_ref()), record.quality());

        let record = reader.next_record().unwrap()?;
        assert_eq!(b"noodles:", record.name());
        assert_eq!(b"2/1", record.description().unwrap());
        assert_eq!(b"TCGA", record.sequence());
        assert_eq!(Some(b"dcba".as_ref()), record.quality());
        assert_eq!(b"noodles:2/1", record.second_description().unwrap());

        let record = reader.next_record().unwrap()?;
        assert_eq!(b"fasta", record.name());
        assert_eq!(
            b"sequence in middle of my fastq",
            record.description().unwrap()
        );
        assert_eq!(b"GATCATGA", record.sequence());
        assert_eq!(None, record.quality());
        assert_eq!(None, record.second_description());

        let record = reader.next_record().unwrap()?;
        assert_eq!(b"1", record.name());
        assert_eq!(None, record.description());
        assert_eq!(b"ACGT", record.sequence());
        assert_eq!(None, record.quality());
        assert_eq!(None, record.second_description());

        let result = reader.next_record();
        println!("{:?}", result);
        assert!(result.is_none());

        Ok(())
    }

    #[test]
    fn bad_record() -> std::io::Result<()> {
        let bad_record = b"";

        let result = Reader::new(&bad_record[..]);
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(std::io::ErrorKind::UnexpectedEof, ""))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
!noodles
ACGT
+
abcd
";

        let mut reader = Reader::new(&bad_record[..])?;
        let result = reader.next_record().unwrap();
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The prefix ('>' or '@') is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
@
ACGT
+
abcd
";

        let mut reader = Reader::new(&bad_record[..])?;
        let result = reader.next_record().unwrap();
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The sequence name is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
@noodles

+
abcd
";

        let mut reader = Reader::new(&bad_record[..])?;
        let result = reader.next_record().unwrap();
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The sequence is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
@noodles
ACTG
abcd
";

        let mut reader = Reader::new(&bad_record[..])?;
        let result = reader.next_record().unwrap();
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The second description is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
@noodles
ACTG
+

";

        let mut reader = Reader::new(&bad_record[..])?;
        let result = reader.next_record().unwrap();
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The quality string is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
>missing sequence

";

        let mut reader = Reader::new(&bad_record[..])?;
        let result = reader.next_record().unwrap();
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The sequence is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
>
ACCTA
";
        let mut reader = Reader::new(&bad_record[..])?;
        let result = reader.next_record().unwrap();
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The sequence name is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
@a
A
";
        let mut reader = Reader::new(&bad_record[..])?;
        let result = reader.next_record().unwrap();
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The last fastq record isn't complete"
                ))
            ),
            format!("{:?}", result.err())
        );

        Ok(())
    }
}
