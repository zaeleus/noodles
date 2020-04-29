use std::io::{self, Write};

use super::Record;

pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: Write,
{
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        writeln!(
            self.inner,
            "{name}\t{len}\t{offset}\t{line_bases}\t{line_width}",
            name = record.name(),
            len = record.len(),
            offset = record.offset(),
            line_bases = record.line_bases(),
            line_width = record.line_width(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());

        let record = Record::new(String::from("sq0"), 10946, 4, 80, 81);
        writer.write_record(&record)?;

        let expected = b"sq0\t10946\t4\t80\t81\n";
        assert_eq!(writer.get_ref(), expected);

        Ok(())
    }
}
