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

        if !intersects_index_record(index_record, self.reference_sequence_id, self.interval) {
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

        let slice = match container.slice(index_record.landmark()) {
            Ok(slice) => slice,
            Err(e) => return Some(Err(e)),
        };

        let (core_data_src, external_data_srcs) = match slice.decode_blocks() {
            Ok(blocks) => blocks,
            Err(e) => return Some(Err(e)),
        };

        let records = slice
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
                        sam::alignment::RecordBuf::try_from_alignment_record(self.header, &record)
                    })
                    .collect::<io::Result<Vec<_>>>()
            });

        let records = match records {
            Ok(records) => records,
            Err(e) => return Some(Err(e)),
        };

        self.records = records.into_iter();

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
                Some(record) => {
                    if intersects(&record, self.reference_sequence_id, self.interval) {
                        return Some(Ok(record));
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

fn intersects(
    record: &sam::alignment::RecordBuf,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> bool {
    match (
        record.reference_sequence_id(),
        record.alignment_start(),
        record.alignment_end(),
    ) {
        (Some(id), Some(start), Some(end)) if id == reference_sequence_id => {
            let alignment_interval = (start..=end).into();
            region_interval.intersects(alignment_interval)
        }
        _ => false,
    }
}

fn intersects_index_record(
    record: &crai::Record,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> bool {
    if record.reference_sequence_id() != Some(reference_sequence_id) {
        return false;
    }

    let Some(start) = record.alignment_start() else {
        return false;
    };

    let Some(end) = start.checked_add(record.alignment_span().saturating_sub(1)) else {
        return true;
    };

    region_interval.intersects((start..=end).into())
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;
    use noodles_sam::alignment::{
        record::cigar::{Op, op::Kind},
        record_buf::Cigar,
    };

    use super::*;

    #[test]
    fn test_intersects() {
        let record = build_record(0, 8, 5);
        let interval =
            Interval::from(Position::try_from(10).unwrap()..=Position::try_from(13).unwrap());

        assert!(intersects(&record, 0, interval));
    }

    #[test]
    fn test_intersects_with_different_reference_sequence() {
        let record = build_record(1, 8, 5);
        let interval =
            Interval::from(Position::try_from(10).unwrap()..=Position::try_from(13).unwrap());

        assert!(!intersects(&record, 0, interval));
    }

    #[test]
    fn test_intersects_index_record() {
        let record = crai::Record::new(Some(0), Position::new(8), 5, 13, 21, 34);
        let interval =
            Interval::from(Position::try_from(10).unwrap()..=Position::try_from(13).unwrap());

        assert!(intersects_index_record(&record, 0, interval));
    }

    #[test]
    fn test_intersects_index_record_with_nonoverlapping_interval() {
        let record = crai::Record::new(Some(0), Position::new(8), 5, 13, 21, 34);
        let interval =
            Interval::from(Position::try_from(13).unwrap()..=Position::try_from(21).unwrap());

        assert!(!intersects_index_record(&record, 0, interval));
    }

    #[test]
    fn test_intersects_index_record_with_different_reference_sequence() {
        let record = crai::Record::new(Some(1), Position::new(8), 5, 13, 21, 34);
        let interval = Interval::from(Position::MIN..=Position::MAX);

        assert!(!intersects_index_record(&record, 0, interval));
    }

    fn build_record(
        reference_sequence_id: usize,
        alignment_start: usize,
        alignment_span: usize,
    ) -> sam::alignment::RecordBuf {
        let cigar: Cigar = [Op::new(Kind::Match, alignment_span)].into_iter().collect();

        sam::alignment::RecordBuf::builder()
            .set_reference_sequence_id(reference_sequence_id)
            .set_alignment_start(Position::try_from(alignment_start).unwrap())
            .set_cigar(cigar)
            .build()
    }
}
