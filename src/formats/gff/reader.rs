use std::io::Read;

pub struct Reader<R: Read> {
    inner: csv::Reader<R>,
}

impl<R: Read> Reader<R> {
    pub fn new(reader: R) -> Reader<R> {
        let inner = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .comment(Some(b'#'))
            .from_reader(reader);

        Reader { inner }
    }

    pub fn records(&mut self) -> csv::StringRecordsIter<R> {
        self.inner.records()
    }
}
