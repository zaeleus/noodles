use std::io::{self, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};

use super::Reader;
use crate::{alignment::Record as _, Header, Record};

pub struct Query<'a, R>
where
    R: Read + Seek,
{
    reader: Reader<csi::io::Query<'a, R>>,
    record: Record,
}

impl<'a, R> Query<'a, R>
where
    R: Read + Seek,
{
    pub(super) fn new(reader: &'a mut bgzf::Reader<R>, chunks: Vec<Chunk>) -> Self {
        Self {
            reader: Reader::new(csi::io::Query::new(reader, chunks)),
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
        match self.next_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

pub struct FilterByRegion<'h, I> {
    records: I,
    header: &'h Header,
    reference_sequence_id: usize,
    interval: Interval,
}

impl<'h, I> FilterByRegion<'h, I>
where
    I: Iterator<Item = io::Result<Record>>,
{
    pub fn new(
        records: I,
        header: &'h Header,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            records,
            header,
            reference_sequence_id,
            interval,
        }
    }
}

impl<'h, I> Iterator for FilterByRegion<'h, I>
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

            match intersects(
                self.header,
                &record,
                self.reference_sequence_id,
                self.interval,
            ) {
                Ok(true) => return Some(Ok(record)),
                Ok(false) => {}
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

fn intersects(
    header: &Header,
    record: &Record,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    match (
        record.reference_sequence_id(header).transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ) {
        (Some(id), Some(start), Some(end)) => {
            let alignment_interval = (start..=end).into();
            Ok(id == reference_sequence_id && region_interval.intersects(alignment_interval))
        }
        _ => Ok(false),
    }
}
