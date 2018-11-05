use std::fs::File;
use std::io::{self, BufReader, Read};
use std::path::Path;

use csv;

use formats::gz::MultiGzDecoder;

pub struct Reader<R: Read> {
    reader: csv::Reader<R>,
}

impl<R: Read> Reader<R> {
    pub fn open<P>(src: P) -> io::Result<Reader<File>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(src)?;
        Ok(Reader::new(file))
    }

    pub fn new(reader: R) -> Reader<R> {
        let reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .comment(Some(b'#'))
            .from_reader(reader);

        Reader { reader }
    }

    pub fn records(&mut self) -> csv::StringRecordsIter<R> {
        self.reader.records()
    }
}

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
