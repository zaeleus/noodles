use std::{
    io::{self, Read},
    vec,
};

use noodles_sam as sam;

use super::{Container, Reader};

/// An iterator over records of a CRAM reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'r, 'h: 'r, R>
where
    R: Read,
{
    reader: &'r mut Reader<R>,
    header: &'h sam::Header,
    container: Container,
    records: vec::IntoIter<sam::alignment::RecordBuf>,
}

impl<'r, 'h: 'r, R> Records<'r, 'h, R>
where
    R: Read,
{
    pub(crate) fn new(reader: &'r mut Reader<R>, header: &'h sam::Header) -> Self {
        Self {
            reader,
            header,
            container: Container::default(),
            records: Vec::new().into_iter(),
        }
    }

    fn read_container_records(&mut self) -> io::Result<bool> {
        if self.reader.read_container(&mut self.container)? == 0 {
            return Ok(true);
        }

        let compression_header = self.container.compression_header()?;

        self.records = self
            .container
            .slices()
            .map(|result| {
                let slice = result?;

                let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

                slice
                    .records(
                        self.reader.reference_sequence_repository.clone(),
                        self.header,
                        &compression_header,
                        &core_data_src,
                        &external_data_srcs,
                    )
                    .and_then(|records| {
                        records
                            .into_iter()
                            .map(|record| {
                                sam::alignment::RecordBuf::try_from_alignment_record(
                                    self.header,
                                    &record,
                                )
                            })
                            .collect::<io::Result<Vec<_>>>()
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

impl<R> Iterator for Records<'_, '_, R>
where
    R: Read,
{
    type Item = io::Result<sam::alignment::RecordBuf>;

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
