mod records;

use std::io;

use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};

use self::records::Records;
use super::Reader;
use crate::{Header, Record, alignment::Record as _};

pub struct Query<'r, 'h: 'r, R> {
    reader: Reader<csi::io::Query<'r, R>>,
    header: &'h Header,
    reference_sequence_id: usize,
    interval: Interval,
}

impl<'r, 'h: 'r, R> Query<'r, 'h, R>
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
        }
    }

    fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        next_record(
            &mut self.reader,
            record,
            self.header,
            self.reference_sequence_id,
            self.interval,
        )
    }

    fn records(self) -> Records<'r, 'h, R> {
        Records::new(self)
    }
}

impl<'r, 'h: 'r, R> IntoIterator for Query<'r, 'h, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    type Item = io::Result<Record>;
    type IntoIter = Records<'r, 'h, R>;

    fn into_iter(self) -> Self::IntoIter {
        self.records()
    }
}

pub(crate) fn intersects(
    header: &Header,
    record: &Record,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    let Some(id) = record.reference_sequence_id(header).transpose()? else {
        return Ok(false);
    };

    if id != reference_sequence_id {
        return Ok(false);
    }

    if interval_is_unbounded(region_interval) {
        Ok(true)
    } else {
        match (
            record.alignment_start().transpose()?,
            record.alignment_end().transpose()?,
        ) {
            (Some(start), Some(end)) => {
                let alignment_interval = (start..=end).into();
                Ok(region_interval.intersects(alignment_interval))
            }
            _ => Ok(false),
        }
    }
}

fn interval_is_unbounded(interval: Interval) -> bool {
    interval.start().is_none() && interval.end().is_none()
}

fn next_record<R>(
    reader: &mut Reader<csi::io::Query<'_, R>>,
    record: &mut Record,
    header: &Header,
    reference_sequence_id: usize,
    interval: Interval,
) -> io::Result<usize>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    loop {
        match reader.read_record(record)? {
            0 => return Ok(0),
            n => {
                if intersects(header, record, reference_sequence_id, interval)? {
                    return Ok(n);
                }
            }
        }
    }
}
