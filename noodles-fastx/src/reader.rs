//! FASTX record

/* std use */

/* crate use */
use crate::record::Record;
use crate::records::Records;

/// Struct to perform reade of fasta or fastq file
pub struct Reader<R>
where
    R: std::io::BufRead,
{
    inner: R,
}

impl<R> Reader<R>
where
    R: std::io::BufRead,
{
    /// Create a new Reader from a [std::io::BufRead]
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Get next record from inner
    pub fn next_record(&mut self) -> Option<std::io::Result<Record>> {
        let mut record = Record::default();

        match record.from_reader(&mut self.inner) {
            Ok(0) => None,
            Ok(_) => Some(Ok(record)),
            Err(e) => Some(Err(e)),
        }
    }

    /// Generate an iterator of Record
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self)
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
";

        let mut reader = Reader::new(&data[..]);
        let mut iterator = reader.records();

        let record = iterator.next().unwrap().unwrap();
        assert_eq!(b"noodles:1/1", record.name());
        assert_eq!(None, record.description());
        assert_eq!(b"AGCT", record.sequence());
        assert_eq!(None, record.second_description());
        assert_eq!(Some(b"abcd".as_ref()), record.quality());

        let record = iterator.next().unwrap().unwrap();
        assert_eq!(b"noodles:", record.name());
        assert_eq!(b"2/1", record.description().unwrap());
        assert_eq!(b"TCGA", record.sequence());
        assert_eq!(Some(b"dcba".as_ref()), record.quality());
        assert_eq!(b"noodles:2/1", record.second_description().unwrap());

        let record = iterator.next().unwrap().unwrap();
        assert_eq!(b"fasta", record.name());
        assert_eq!(
            b"sequence in middle of my fastq",
            record.description().unwrap()
        );
        assert_eq!(b"GATCATGA", record.sequence());
        assert_eq!(None, record.quality());
        assert_eq!(None, record.second_description());

        assert!(iterator.next().is_none());

        Ok(())
    }

    #[test]
    fn bad_record() -> std::io::Result<()> {
        let bad_record = b"\
!noodles
ACGT
+
abcd
";

        let mut reader = Reader::new(&bad_record[..]);
        let mut iterator = reader.records();

        let result = iterator.next().unwrap();
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

        Ok(())
    }
}
