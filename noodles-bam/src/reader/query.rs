use std::io::{self, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_core::{region::Interval, Position};
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_sam::{self as sam, alignment::Record as _};

use super::Reader;
use crate::Record;

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
    record: Record,
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
            record: Record::default(),
        }
    }

    fn next_record(&mut self) -> io::Result<Option<Record>> {
        self.reader.read_record(&mut self.record).map(|n| match n {
            0 => None,
            _ => Some(self.record.clone()),
        })
    }
}

impl<'a, R> Iterator for Query<'a, R>
where
    R: Read + Seek,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.next_record() {
                Ok(Some(record)) => {
                    match intersects(
                        &record,
                        self.header,
                        self.reference_sequence_id,
                        self.interval,
                    ) {
                        Ok(true) => return Some(Ok(record)),
                        Ok(false) => {}
                        Err(e) => return Some(Err(e)),
                    }
                }
                Ok(None) => return None,
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

fn intersects(
    record: &Record,
    header: &sam::Header,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    match (
        record
            .reference_sequence_id()
            .map(usize::try_from)
            .transpose()?,
        record
            .alignment_start()
            .map(Position::try_from)
            .transpose()?,
        record.alignment_end(header).transpose()?,
    ) {
        (Some(id), Some(start), Some(end)) => {
            let alignment_interval = (start..=end).into();
            Ok(id == reference_sequence_id && region_interval.intersects(alignment_interval))
        }
        _ => Ok(false),
    }
}
