use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use crate::formats::fastq::Record;
use crate::formats::gz::MultiGzDecoder;

pub struct Reader<R: BufRead> {
    reader: R,
}

impl<R: BufRead> Reader<R> {
    pub fn new(reader: R) -> Reader<R> {
        Reader { reader }
    }

    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        record.clear();

        let mut len = read_line(&mut self.reader, record.name_mut())?;

        if len > 0 {
            len += read_line(&mut self.reader, record.sequence_mut())?;
            len += read_line(&mut self.reader, record.plus_line_mut())?;
            len += read_line(&mut self.reader, record.quality_mut())?;
        }

        Ok(len)
    }
}

pub fn open<P>(src: P) -> io::Result<Reader<Box<dyn BufRead>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let extenstion = path.extension();
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    match extenstion.and_then(|ext| ext.to_str()) {
        Some("gz") => {
            let decoder = MultiGzDecoder::new(reader);
            Ok(Reader::new(Box::new(BufReader::new(decoder))))
        },
        _ => {
            Ok(Reader::new(Box::new(reader)))
        }
    }
}

fn read_line<R: BufRead>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize> {
    let result = reader.read_until(b'\n', buf);

    // Chomp newline.
    if result.is_ok() {
        let len = buf.len();

        if len > 0 {
            buf.truncate(len - 1);
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;

    use crate::formats::fastq::Record;
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
        assert_eq!(record, Record::new("@noodles:1/1", "AGCT", "+", "abcd"));

        let len = reader.read_record(&mut record).unwrap();
        assert_eq!(len, 25);
        assert_eq!(record, Record::new("@noodles:2/1", "TCGA", "+", "dcba"));

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
