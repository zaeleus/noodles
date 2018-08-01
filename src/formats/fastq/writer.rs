use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use formats::bam::Record;
use formats::bam::sequence::Complement;

pub struct Writer<W> {
    writer: W,
}

impl<W: Write> Writer<W> {
    pub fn create<P>(dst: P) -> io::Result<Writer<BufWriter<File>>> where P: AsRef<Path> {
        let file = File::create(dst)?;
        let writer = BufWriter::new(file);
        Ok(Writer::new(writer))
    }

    pub fn new(writer: W) -> Writer<W> {
        Writer { writer }
    }

    pub fn write_bam_record(&mut self, record: &Record) -> io::Result<()> {
        let name = record.read_name();

        let (seq, qual): (String, String) = if record.flag().is_reverse() {
            (Complement::new(record.sequence().symbols().rev()).collect(),
                record.quality().chars().rev().collect())
        } else {
            (record.sequence().symbols().collect(),
                record.quality().chars().collect())
        };

        self.writer.write(b"@")?;
        self.writer.write(name.as_bytes())?;
        self.writer.write(b"\n")?;
        self.writer.write(seq.as_bytes())?;
        self.writer.write(b"\n")?;
        self.writer.write(b"+\n")?;
        self.writer.write(qual.as_bytes())?;
        self.writer.write(b"\n")?;

        Ok(())
    }
}
