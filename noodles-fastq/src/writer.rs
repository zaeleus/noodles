use std::io::{self, Write};

use super::Record;

pub struct Writer<W> {
    writer: W,
}

impl<W: Write> Writer<W> {
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        self.writer.write_all(record.name())?;
        self.writer.write_all(b"\n")?;
        self.writer.write_all(record.sequence())?;
        self.writer.write_all(b"\n")?;
        self.writer.write_all(record.plus_line())?;
        self.writer.write_all(b"\n")?;
        self.writer.write_all(record.quality())?;
        self.writer.write_all(b"\n")?;

        Ok(())
    }
}
