use std::io::Read;

pub struct Reader<R> {
    reader: csv::Reader<R>,
}

impl<R> Reader<R>
where
    R: Read,
{
    pub fn new(inner: R) -> Self {
        let reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(inner);

        Self { reader }
    }

    pub fn records(&mut self) -> csv::StringRecordsIter<R> {
        self.reader.records()
    }
}
