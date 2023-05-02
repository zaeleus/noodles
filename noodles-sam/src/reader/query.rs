use std::io::{self, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::{self as csi, index::reference_sequence::bin::Chunk};

use super::Reader;
use crate::{alignment::Record, Header};

pub struct Query<'a, R>
where
    R: Read + Seek,
{
    reader: Reader<csi::io::Query<'a, R>>,
    header: &'a Header,
    record: Record,
}

impl<'a, R> Query<'a, R>
where
    R: Read + Seek,
{
    pub(super) fn new(
        reader: &'a mut bgzf::Reader<R>,
        header: &'a Header,
        chunks: Vec<Chunk>,
    ) -> Self {
        Self {
            reader: Reader::new(csi::io::Query::new(reader, chunks)),
            header,
            record: Record::default(),
        }
    }

    fn next_record(&mut self) -> io::Result<Option<Record>> {
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
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

pub struct FilterByRegion<I> {
    records: I,
    reference_sequence_id: usize,
    interval: Interval,
}

impl<I> FilterByRegion<I>
where
    I: Iterator<Item = io::Result<Record>>,
{
    pub fn new(records: I, reference_sequence_id: usize, interval: Interval) -> Self {
        Self {
            records,
            reference_sequence_id,
            interval,
        }
    }
}

impl<I> Iterator for FilterByRegion<I>
where
    I: Iterator<Item = io::Result<Record>>,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let record = match self.records.next()? {
                Ok(record) => record,
                Err(e) => return Some(Err(e)),
            };

            if let (Some(id), Some(start), Some(end)) = (
                record.reference_sequence_id(),
                record.alignment_start(),
                record.alignment_end(),
            ) {
                let alignment_interval = Interval::from(start..=end);

                if id == self.reference_sequence_id && self.interval.intersects(alignment_interval)
                {
                    return Some(Ok(record));
                }
            }
        }
    }
}
