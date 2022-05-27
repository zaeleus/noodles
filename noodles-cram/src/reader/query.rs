use std::{
    io::{self, Read, Seek, SeekFrom},
    ops::RangeBounds,
    slice, vec,
};

use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::Reader;
use crate::{crai, Record};

/// An iterator over records that intersect a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'a, R, B>
where
    R: Read + Seek,
    B: RangeBounds<Position> + Copy,
{
    reader: &'a mut Reader<R>,

    reference_sequence_repository: &'a fasta::Repository,
    header: &'a sam::Header,

    index: slice::Iter<'a, crai::Record>,

    reference_sequence_id: usize,
    interval: B,

    records: vec::IntoIter<Record>,
}

impl<'a, R, B> Query<'a, R, B>
where
    R: Read + Seek,
    B: RangeBounds<Position> + Copy,
{
    pub(super) fn new(
        reader: &'a mut Reader<R>,
        reference_sequence_repository: &'a fasta::Repository,
        header: &'a sam::Header,
        index: &'a crai::Index,
        reference_sequence_id: usize,
        interval: B,
    ) -> Self {
        Self {
            reader,

            reference_sequence_repository,
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

        let container = match self.reader.read_data_container() {
            Ok(Some(c)) => c,
            Ok(None) => return None,
            Err(e) => return Some(Err(e)),
        };

        let records = container
            .slices()
            .iter()
            .map(|slice| {
                let compression_header = container.compression_header();

                slice.records(compression_header).and_then(|mut records| {
                    slice.resolve_records(
                        self.reference_sequence_repository,
                        self.header,
                        compression_header,
                        &mut records,
                    )?;

                    Ok(records)
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

impl<'a, R, B> Iterator for Query<'a, R, B>
where
    R: Read + Seek,
    B: RangeBounds<Position> + Copy,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.records.next() {
                Some(r) => {
                    if let (Some(start), Some(end)) = (r.alignment_start(), r.alignment_end()) {
                        if in_interval(start, end, self.interval) {
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

fn in_interval<B>(alignment_start: Position, alignment_end: Position, region_interval: B) -> bool
where
    B: RangeBounds<Position>,
{
    use std::ops::Bound;

    let a = match region_interval.start_bound() {
        Bound::Included(start) => *start <= alignment_end,
        Bound::Excluded(start) => *start < alignment_end,
        Bound::Unbounded => true,
    };

    let b = match region_interval.end_bound() {
        Bound::Included(end) => alignment_start <= *end,
        Bound::Excluded(end) => alignment_start < *end,
        Bound::Unbounded => true,
    };

    a && b
}
