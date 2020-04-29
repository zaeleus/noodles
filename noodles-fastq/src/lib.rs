pub use self::reader::Reader;
pub use self::record::Record;
pub use self::writer::Writer;

pub mod reader;
pub mod record;
pub mod writer;

use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::Path,
};

use flate2::{bufread::MultiGzDecoder, write::GzEncoder, Compression};

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
        }
        _ => Ok(Reader::new(Box::new(reader))),
    }
}
