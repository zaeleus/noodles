use std::{
    io::{self, Read},
    vec,
};

use noodles_fasta as fasta;

use super::Reader;
use crate::Record;

/// An iterator over records of a CRAM reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'a, 'b, R>
where
    R: Read,
{
    reader: &'a mut Reader<R>,
    reference_sequences: &'b [fasta::Record],
    records: vec::IntoIter<Record>,
}

impl<'a, 'b, R> Records<'a, 'b, R>
where
    R: Read,
{
    pub(crate) fn new(reader: &'a mut Reader<R>, reference_sequences: &'b [fasta::Record]) -> Self {
        Self {
            reader,
            reference_sequences,
            records: Vec::new().into_iter(),
        }
    }

    fn read_container_records(&mut self) -> io::Result<bool> {
        let container = match self.reader.read_data_container()? {
            Some(c) => c,
            None => return Ok(true),
        };

        self.records = container
            .slices()
            .iter()
            .map(|slice| {
                let compression_header = container.compression_header();

                slice.records(compression_header).and_then(|records| {
                    slice.resolve_records(self.reference_sequences, compression_header, records)
                })
            })
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flatten()
            .collect::<Vec<_>>()
            .into_iter();

        Ok(false)
    }
}

impl<'a, 'b, R> Iterator for Records<'a, 'b, R>
where
    R: Read,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.records.next() {
                Some(r) => return Some(Ok(r)),
                None => match self.read_container_records() {
                    Ok(true) => return None,
                    Ok(false) => {}
                    Err(e) => return Some(Err(e)),
                },
            }
        }
    }
}
