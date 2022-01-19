use std::{
    io::{self, Read, Seek},
    ops::{Bound, RangeBounds},
};

use noodles_bgzf as bgzf;
use noodles_csi::index::reference_sequence::bin::Chunk;

use crate::Record;

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

    chunks: Vec<Chunk>,
    i: usize,

    chromosome_id: usize,
    start: i32,
    end: i32,

    state: State,
    record: Record,
}

impl<'a, R> Query<'a, R>
where
    R: Read + Seek,
{
    pub(crate) fn new<B>(
        reader: &'a mut Reader<bgzf::Reader<R>>,
        chunks: Vec<Chunk>,
        chromosome_id: usize,
        interval: B,
    ) -> Self
    where
        B: RangeBounds<i32>,
    {
        let (start, end) = resolve_interval(interval);

        Self {
            reader,

            chunks,
            i: 0,

            chromosome_id,
            start,
            end,

            state: State::Seek,
            record: Record::default(),
        }
    }

    fn read_record(&mut self) -> io::Result<Option<Record>> {
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
            match self.state {
                State::Seek => {
                    self.state = match next_chunk(&self.chunks, &mut self.i) {
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

                        match intersects(&record, self.chromosome_id, self.start, self.end) {
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

fn resolve_interval<B>(interval: B) -> (i32, i32)
where
    B: RangeBounds<i32>,
{
    match (interval.start_bound(), interval.end_bound()) {
        (Bound::Included(s), Bound::Included(e)) => (*s, *e),
        (Bound::Included(s), Bound::Unbounded) => (*s, i32::MAX),
        (Bound::Unbounded, Bound::Unbounded) => (1, i32::MAX),
        _ => todo!(),
    }
}

fn next_chunk(chunks: &[Chunk], i: &mut usize) -> Option<Chunk> {
    let chunk = chunks.get(*i).copied();
    *i += 1;
    chunk
}

fn intersects(
    record: &Record,
    chromosome_id: usize,
    interval_start: i32,
    interval_end: i32,
) -> io::Result<bool> {
    let id = record.chromosome_id();

    let start = i32::from(record.position());
    let end = record.end().map(i32::from)?;

    Ok(id == chromosome_id && in_interval(start, end, interval_start, interval_end))
}

fn in_interval(a_start: i32, a_end: i32, b_start: i32, b_end: i32) -> bool {
    a_start <= b_end && b_start <= a_end
}
