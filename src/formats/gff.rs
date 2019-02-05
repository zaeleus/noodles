pub use self::{attributes::Attributes, reader::Reader, record::Record, strand::Strand};

pub mod attributes;
pub mod reader;
pub mod record;
pub mod strand;

use std::{
    fs::File,
    io::{self, BufReader, Read},
    path::Path,
};

use crate::formats::gz::MultiGzDecoder;

pub fn open<P>(src: P) -> io::Result<Reader<Box<dyn Read>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let extension = path.extension();
    let file = File::open(path)?;

    match extension.and_then(|ext| ext.to_str()) {
        Some("gz") => {
            let reader = BufReader::new(file);
            let decoder = MultiGzDecoder::new(reader);
            Ok(Reader::new(Box::new(decoder)))
        },
        _ => {
            Ok(Reader::new(Box::new(file)))
        }
    }
}
