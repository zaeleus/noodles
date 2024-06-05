use std::io;

use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};

use super::Reader;
use crate::{alignment::Record as _, Header, Record};

pub struct Query<'r, 'h, R> {
    reader: Reader<csi::io::Query<'r, R>>,
    header: &'h Header,
    reference_sequence_id: usize,
    interval: Interval,
    record: Record,
}

impl<'r, 'h, R> Query<'r, 'h, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    pub(super) fn new(
        reader: &'r mut R,
        chunks: Vec<Chunk>,
        header: &'h Header,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader: Reader::new(csi::io::Query::new(reader, chunks)),
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

impl<'r, 'h, R> Iterator for Query<'r, 'h, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let record = match self.next_record() {
                Ok(Some(record)) => record,
                Ok(None) => return None,
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
