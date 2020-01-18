use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use flate2::{write::GzEncoder, Compression};

use super::Record;

pub struct Writer<W> {
    writer: W,
}

impl<W: Write> Writer<W> {
    pub fn create<P>(dst: P) -> io::Result<Writer<BufWriter<File>>>
    where
        P: AsRef<Path>,
    {
        let file = File::create(dst)?;
        let writer = BufWriter::new(file);
        Ok(Writer::new(writer))
    }

    pub fn new(writer: W) -> Writer<W> {
        Writer { writer }
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
        }
        _ => Ok(Writer::new(Box::new(writer))),
    }
}
