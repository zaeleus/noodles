use std::io;

use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_sam::alignment::Record as _;
use noodles_sam::{self as sam, alignment::RecordBuf};

use super::{Reader, RecordBufs, Records};
use crate::Record;

/// An iterator over records of a BAM reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'a, R> {
    reader: Reader<csi::io::Query<'a, R>>,
    reference_sequence_id: usize,
    interval: Interval,
    record: Record,
}

impl<'a, R> Query<'a, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    pub(super) fn new(
        reader: &'a mut R,
        chunks: Vec<Chunk>,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader: Reader::from(csi::io::Query::new(reader, chunks)),
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

    /// Create a `QuriedReader`, which behaves like other `Reader`s.
    pub fn into_reader(self) -> QueriedReader<'a, R> {
        QueriedReader {
            reader: self.reader,
            reference_sequence_id: self.reference_sequence_id,
            interval: self.interval,
        }
    }

}

impl<'a, R> Iterator for Query<'a, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.next_record() {
                Ok(Some(record)) => {
                    match intersects(&record, self.reference_sequence_id, self.interval) {
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

pub(crate) fn intersects(
    record: &Record,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    match (
        record.reference_sequence_id().transpose()?,
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

/// A queried BAM reader.
pub struct QueriedReader<'a, R> {
    reader: Reader<csi::io::Query<'a, R>>,
    reference_sequence_id: usize,
    interval: Interval,
    // record: Record,
}

impl<'a, R> QueriedReader<'a, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    /// Reads a record into an alignment record buffer.
    pub fn read_record_buf(
        &mut self,
        header: &sam::Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        loop {
            match self.reader.read_record_buf(header, record) {
                Ok(0) => return Ok(0),
                Ok(n) => match intersects_buf(&record, self.reference_sequence_id, self.interval) {
                    true => return Ok(n),
                    false => {}
                },
                Err(e) => return Err(e),
            }
        }
    }

    /// Reads a record.
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        loop {
            match self.reader.read_record(record) {
                Ok(0) => return Ok(0),
                Ok(n) => match intersects(&record, self.reference_sequence_id, self.interval) {
                    Ok(true) => return Ok(n),
                    Ok(false) => {}
                    Err(e) => return Err(e),
                }
                Err(e) => return Err(e),
            }
        }
    }

    /// Returns an iterator over alignment record buffers starting from the current stream
    /// position.
    pub fn record_bufs(
        &'a mut self,
        header: &'a sam::Header,
    ) -> RecordBufs<'_, csi::io::Query<'a, R>> {
        self.reader.record_bufs(header)
    }

    /// Returns an iterator over records.
    pub fn records(&mut self) -> Records<'_, csi::io::Query<'a, R>> {
        self.reader.records()
    }
}

fn intersects_buf(
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