use std::io::{self, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_sam::{self as sam, alignment::RecordBuf};

use super::Reader;

/// An iterator over records of a BAM reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'a, R>
where
    R: Read + Seek,
{
    reader: Reader<csi::io::Query<'a, R>>,
    header: &'a sam::Header,
    reference_sequence_id: usize,
    interval: Interval,
    record: RecordBuf,
}

impl<'a, R> Query<'a, R>
where
    R: Read + Seek,
{
    pub(super) fn new(
        reader: &'a mut bgzf::Reader<R>,
        header: &'a sam::Header,
        chunks: Vec<Chunk>,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader: Reader::from(csi::io::Query::new(reader, chunks)),
            header,
            reference_sequence_id,
            interval,
            record: RecordBuf::default(),
        }
    }

    fn next_record(&mut self) -> io::Result<Option<RecordBuf>> {
        self.reader
            .read_record(self.header, &mut self.record)
            .map(|n| match n {
                0 => None,
                _ => Some(self.record.clone()),
            })
    }
}

impl<'a, R> Iterator for Query<'a, R>
where
    R: Read + Seek,
{
    type Item = io::Result<RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.next_record() {
                Ok(Some(record)) => {
                    if intersects(&record, self.reference_sequence_id, self.interval) {
                        return Some(Ok(record));
                    }
                }
                Ok(None) => return None,
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

pub(crate) fn intersects(
    record: &RecordBuf,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> bool {
    match (
        record.reference_sequence_id(),
        record.alignment_start(),
        record.alignment_end(),
    ) {
        (Some(id), Some(start), Some(end)) => {
            let alignment_interval = (start..=end).into();
            id == reference_sequence_id && region_interval.intersects(alignment_interval)
        }
        _ => false,
    }
}
