use std::{
    io::{self, Read, Seek},
    vec,
};

use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::index::reference_sequence::bin::Chunk;

use super::Reader;
use crate::{alignment::Record, Header};

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

pub struct Query<'a, R>
where
    R: Read + Seek,
{
    reader: &'a mut Reader<bgzf::Reader<R>>,
    header: &'a Header,
    chunks: vec::IntoIter<Chunk>,
    state: State,
    record: Record,
}

impl<'a, R> Query<'a, R>
where
    R: Read + Seek,
{
    pub(super) fn new(
        reader: &'a mut Reader<bgzf::Reader<R>>,
        header: &'a Header,
        chunks: Vec<Chunk>,
    ) -> Self {
        Self {
            reader,
            header,
            chunks: chunks.into_iter(),
            state: State::Seek,
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
        loop {
            match self.state {
                State::Seek => {
                    self.state = match self.chunks.next() {
                        Some(chunk) => {
                            if let Err(e) = self.reader.seek(chunk.start()) {
                                return Some(Err(e));
                            }

                            State::Read(chunk.end())
                        }
                        None => State::Done,
                    }
                }
                State::Read(chunk_end) => match self.next_record() {
                    Ok(Some(record)) => {
                        if self.reader.get_ref().virtual_position() >= chunk_end {
                            self.state = State::Seek;
                        }

                        return Some(Ok(record));
                    }
                    Ok(None) => self.state = State::Seek,
                    Err(e) => return Some(Err(e)),
                },
                State::Done => return None,
            }
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
