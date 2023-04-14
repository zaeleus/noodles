use std::{
    io::{self, Read, Seek},
    vec,
};

use noodles_bgzf as bgzf;
use noodles_core::{region::Interval, Position};
use noodles_csi::index::reference_sequence::bin::Chunk;

use crate::lazy;

use super::Reader;

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

/// An iterator over records of a BCF reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'a, R>
where
    R: Read + Seek,
{
    reader: &'a mut Reader<bgzf::Reader<R>>,

    chunks: vec::IntoIter<Chunk>,

    chromosome_id: usize,
    interval: Interval,

    state: State,
    record: lazy::Record,
}

impl<'a, R> Query<'a, R>
where
    R: Read + Seek,
{
    pub(crate) fn new(
        reader: &'a mut Reader<bgzf::Reader<R>>,
        chunks: Vec<Chunk>,
        chromosome_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader,

            chunks: chunks.into_iter(),

            chromosome_id,
            interval,

            state: State::Seek,
            record: lazy::Record::default(),
        }
    }

    fn read_record(&mut self) -> io::Result<Option<lazy::Record>> {
        self.reader
            .read_lazy_record(&mut self.record)
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
    type Item = io::Result<lazy::Record>;

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
                State::Read(chunk_end) => match self.read_record() {
                    Ok(Some(record)) => {
                        if self.reader.virtual_position() >= chunk_end {
                            self.state = State::Seek;
                        }

                        match intersects(&record, self.chromosome_id, self.interval) {
                            Ok(true) => return Some(Ok(record)),
                            Ok(false) => {}
                            Err(e) => return Some(Err(e)),
                        }
                    }
                    Ok(None) => self.state = State::Seek,
                    Err(e) => return Some(Err(e)),
                },
                State::Done => return None,
            }
        }
    }
}

pub(crate) fn intersects(
    record: &lazy::Record,
    chromosome_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    let id = record.chromosome_id();

    let start = Position::try_from(usize::from(record.position()))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let end = record.end().map(usize::from).and_then(|n| {
        Position::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let record_interval = Interval::from(start..=end);

    Ok(id == chromosome_id && record_interval.intersects(region_interval))
}
