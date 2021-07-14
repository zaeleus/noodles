use std::{
    convert::TryFrom,
    io::{self, Read},
    vec,
};

use crate::{DataContainer, Record};

use super::Reader;

/// An iterator over records of a CRAM reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'a, R>
where
    R: Read,
{
    reader: &'a mut Reader<R>,
    records: vec::IntoIter<Record>,
}

impl<'a, R> Records<'a, R>
where
    R: Read,
{
    pub(crate) fn new(reader: &'a mut Reader<R>) -> Records<'_, R> {
        Self {
            reader,
            records: Vec::new().into_iter(),
        }
    }

    fn read_container_records(&mut self) -> io::Result<bool> {
        let container = self.reader.read_container()?;

        if container.is_eof() {
            return Ok(true);
        }

        let data_container = DataContainer::try_from(container)?;

        self.records = data_container
            .slices()
            .iter()
            .map(|slice| {
                slice
                    .records(data_container.compression_header())
                    .map(|r| slice.resolve_mates(r))
            })
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flatten()
            .collect::<Vec<_>>()
            .into_iter();

        Ok(false)
    }
}

impl<'a, R> Iterator for Records<'a, R>
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
