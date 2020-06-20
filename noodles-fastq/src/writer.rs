use std::io::{self, Write};

use super::Record;

pub struct Writer<W> {
    inner: W,
}

impl<W: Write> Writer<W> {
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        self.inner.write_all(record.read_name())?;
        self.inner.write_all(b"\n")?;
        self.inner.write_all(record.sequence())?;
        self.inner.write_all(b"\n+\n")?;
        self.inner.write_all(record.quality_scores())?;
        self.inner.write_all(b"\n")?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());

        let record = Record::new("r0", "ATCG", "NDLS");
        writer.write_record(&record)?;

        let expected = b"r0\nATCG\n+\nNDLS\n";
        assert_eq!(writer.get_ref(), expected);

        Ok(())
    }
}
