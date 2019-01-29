use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use flate2::Compression;
use flate2::write::GzEncoder;

use crate::formats::{bam, fastq};
use crate::formats::bam::sequence::Complement;

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

    pub fn write_bam_record(&mut self, record: &bam::Record) -> io::Result<()> {
        let name = record.read_name();

        let (seq, qual): (String, String) = if record.flag().is_reverse() {
            (Complement::new(record.sequence().symbols().rev()).collect(),
                record.quality().chars().rev().collect())
        } else {
            (record.sequence().symbols().collect(),
                record.quality().chars().collect())
        };

        self.writer.write_all(b"@")?;
        self.writer.write_all(name)?;
        self.writer.write_all(b"\n")?;
        self.writer.write_all(seq.as_bytes())?;
        self.writer.write_all(b"\n")?;
        self.writer.write_all(b"+\n")?;
        self.writer.write_all(qual.as_bytes())?;
        self.writer.write_all(b"\n")?;

        Ok(())
    }

    pub fn write_record(&mut self, record: &fastq::Record) -> io::Result<()> {
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

pub fn create<P>(dst: P) -> io::Result<Writer<Box<dyn Write>>>
where
    P: AsRef<Path>,
{
    let path = dst.as_ref();
    let extension = path.extension();
    let file = File::create(path)?;
    let writer = BufWriter::new(file);

    match extension.and_then(|ext| ext.to_str()) {
        Some("gz") => {
            let level = Compression::default();
            let encoder = GzEncoder::new(writer, level);
            Ok(Writer::new(Box::new(encoder)))
        },
        _ => {
            Ok(Writer::new(Box::new(writer)))
        },
    }
}
