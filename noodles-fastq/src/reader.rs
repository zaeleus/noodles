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
        record.clear();

        let mut len = read_line(&mut self.inner, record.read_name_mut())?;

        if len > 0 {
            len += read_line(&mut self.inner, record.sequence_mut())?;
            len += read_line(&mut self.inner, &mut Vec::new())?;
            len += read_line(&mut self.inner, record.quality_scores_mut())?;
        }

        Ok(len)
    }
}

fn read_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    let result = reader.read_until(b'\n', buf);
    buf.pop();
    result
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;

    use crate::Record;

    use super::{read_line, Reader};

    #[test]
    fn test_read_record() {
        let data = "\
@noodles:1/1
AGCT
+
abcd
@noodles:2/1
TCGA
+
dcba
";

        let mut reader = Reader::new(data.as_bytes());
        let mut record = Record::default();

        let len = reader.read_record(&mut record).unwrap();
        assert_eq!(len, 25);
        assert_eq!(record, Record::new("@noodles:1/1", "AGCT", "abcd"));

        let len = reader.read_record(&mut record).unwrap();
        assert_eq!(len, 25);
        assert_eq!(record, Record::new("@noodles:2/1", "TCGA", "dcba"));

        let len = reader.read_record(&mut record).unwrap();
        assert_eq!(len, 0);
    }

    #[test]
    fn test_read_line() {
        let data = "@fqlib\nAGCT\n";
        let mut reader = BufReader::new(data.as_bytes());

        let mut buf = Vec::new();
        let len = read_line(&mut reader, &mut buf).unwrap();
        assert_eq!(len, 7);
        assert_eq!(buf, b"@fqlib");

        buf.clear();
        let len = read_line(&mut reader, &mut buf).unwrap();
        assert_eq!(len, 5);
        assert_eq!(buf, b"AGCT");

        buf.clear();
        let len = read_line(&mut reader, &mut buf).unwrap();
        assert_eq!(len, 0);
    }
}
