use std::{
    io::{self, Read, Seek, SeekFrom},
    slice, vec,
};

use noodles_core::region::Interval;
use noodles_sam as sam;

use super::{Container, Reader};
use crate::crai;

/// An iterator over records that intersect a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'r, 'h: 'r, 'i: 'r, R>
where
    R: Read + Seek,
{
    reader: &'r mut Reader<R>,

    header: &'h sam::Header,

    index: slice::Iter<'i, crai::Record>,

    reference_sequence_id: usize,
    interval: Interval,

    records: vec::IntoIter<sam::alignment::RecordBuf>,
}

impl<'r, 'h: 'r, 'i: 'r, R> Query<'r, 'h, 'i, R>
where
    R: Read + Seek,
{
    pub(super) fn new(
        reader: &'r mut Reader<R>,
        header: &'h sam::Header,
        index: &'i crai::Index,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader,

            header,

            index: index.iter(),

            reference_sequence_id,
            interval,

            records: Vec::new().into_iter(),
        }
    }

    fn read_next_container(&mut self) -> Option<io::Result<()>> {
        let index_record = self.index.next()?;

        if index_record.reference_sequence_id() != Some(self.reference_sequence_id) {
            return Some(Ok(()));
        }

        if let Err(e) = self.reader.seek(SeekFrom::Start(index_record.offset())) {
            return Some(Err(e));
        }

        let mut container = Container::default();

        match self.reader.read_container(&mut container) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(e)),
        };

        let compression_header = match container.compression_header() {
            Ok(compression_header) => compression_header,
            Err(e) => return Some(Err(e)),
        };

        let records = container
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
            .collect::<Result<Vec<_>, _>>();

        let records = match records {
            Ok(records) => records,
            Err(e) => return Some(Err(e)),
        };

        self.records = records
            .into_iter()
            .flatten()
            .collect::<Vec<_>>()
            .into_iter();

        Some(Ok(()))
    }
}

impl<R> Iterator for Query<'_, '_, '_, R>
where
    R: Read + Seek,
{
    type Item = io::Result<sam::alignment::RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.records.next() {
                Some(r) => {
                    if let (Some(start), Some(end)) = (r.alignment_start(), r.alignment_end()) {
                        let alignment_interval = (start..=end).into();

                        if self.interval.intersects(alignment_interval) {
                            return Some(Ok(r));
                        }
                    }
                }
                None => match self.read_next_container() {
                    Some(Ok(())) => {}
                    Some(Err(e)) => return Some(Err(e)),
                    None => return None,
                },
            }
        }
    }
}
