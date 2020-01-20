use std::io::Read;

pub struct Reader<R: Read> {
    inner: csv::Reader<R>,
}

impl<R: Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        let inner = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .comment(Some(b'#'))
            .from_reader(reader);

        Self { inner }
    }

    pub fn records(&mut self) -> csv::StringRecordsIter<R> {
        self.inner.records()
    }
}
