use std::io::{self, Read};

use noodles_fasta as fasta;
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
}

impl<'r, 'h: 'r, R> Records<'r, 'h, R>
where
    R: Read,
{
    pub(crate) fn new(reader: &'r mut Reader<R>, header: &'h sam::Header) -> Self {
        Self { reader, header }
    }
}

impl<R> Iterator for Records<'_, '_, R>
where
    R: Read,
{
    type Item = io::Result<sam::alignment::RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = sam::alignment::RecordBuf::default();

        match self.reader.read_record_buf(self.header, &mut record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(record)),
            Err(e) => Some(Err(e)),
        }
    }
}

/// Decodes all records of a container into owned alignment records.
pub(super) fn decode_container_records(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    container: &Container,
) -> io::Result<Vec<sam::alignment::RecordBuf>> {
    let compression_header = container.compression_header()?;

    let records = container
        .slices()
        .map(|result| {
            let slice = result?;

            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

            slice
                .records(
                    reference_sequence_repository.clone(),
                    header,
                    &compression_header,
                    &core_data_src,
                    &external_data_srcs,
                )
                .and_then(|records| {
                    records
                        .into_iter()
                        .map(|record| {
                            sam::alignment::RecordBuf::try_from_alignment_record(header, &record)
                        })
                        .collect::<io::Result<Vec<_>>>()
                })
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(records.into_iter().flatten().collect())
}
