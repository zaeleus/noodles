use std::{
    io::{self, Read, Seek},
    ops::{Bound, RangeBounds},
};

use noodles_bgzf::{self as bgzf, VirtualPosition};
use noodles_core::Position;
use noodles_csi::index::reference_sequence::bin::Chunk;
use noodles_sam::AlignmentRecord;

use crate::Record;

use super::Reader;

enum State {
    Seek,
    Read(VirtualPosition),
    Done,
}

/// An iterator over records of a BAM reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'a, R, B>
where
    R: Read + Seek,
    B: RangeBounds<Position> + Copy,
{
    reader: &'a mut Reader<bgzf::Reader<R>>,

    chunks: Vec<Chunk>,
    i: usize,

    reference_sequence_id: usize,
    interval: B,

    state: State,
    record: Record,
}

impl<'a, R, B> Query<'a, R, B>
where
    R: Read + Seek,
    B: RangeBounds<Position> + Copy,
{
    pub(super) fn new(
        reader: &'a mut Reader<bgzf::Reader<R>>,
        chunks: Vec<Chunk>,
        reference_sequence_id: usize,
        interval: B,
    ) -> Self {
        Self {
            reader,

            chunks,
            i: 0,

            reference_sequence_id,
            interval,

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

impl<'a, R, B> Iterator for Query<'a, R, B>
where
    R: Read + Seek,
    B: RangeBounds<Position> + Copy,
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

                        match intersects(&record, self.reference_sequence_id, self.interval) {
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

pub(crate) fn next_chunk(chunks: &[Chunk], i: &mut usize) -> Option<Chunk> {
    let chunk = chunks.get(*i).copied();
    *i += 1;
    chunk
}

pub(crate) fn intersects<B>(
    record: &Record,
    reference_sequence_id: usize,
    interval: B,
) -> io::Result<bool>
where
    B: RangeBounds<Position>,
{
    let id = match record.reference_sequence_id() {
        Some(id) => id,
        None => return Ok(false),
    };

    let start = record.alignment_start().expect("missing alignment start");
    let end = record.alignment_end().expect("missing alignment end");

    Ok(id == reference_sequence_id && in_interval(start, end, interval))
}

fn in_interval<B>(alignment_start: Position, alignment_end: Position, region_interval: B) -> bool
where
    B: RangeBounds<Position>,
{
    let a = match region_interval.start_bound() {
        Bound::Included(start) => *start <= alignment_end,
        Bound::Excluded(start) => *start < alignment_end,
        Bound::Unbounded => true,
    };

    let b = match region_interval.end_bound() {
        Bound::Included(end) => alignment_start <= *end,
        Bound::Excluded(end) => alignment_start < *end,
        Bound::Unbounded => true,
    };

    a && b
}
