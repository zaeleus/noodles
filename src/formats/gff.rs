use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

use csv;

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
