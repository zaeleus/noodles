use std::{
    io::{self, Read, Seek},
    ops::{Bound, RangeBounds},
};

use noodles_bgzf::{self as bgzf, VirtualPosition};
use noodles_core::{region::Interval, Position};
use noodles_csi::index::reference_sequence::bin::Chunk;
use noodles_sam::alignment::Record;

use super::Reader;

enum State {
    Seek,
    Read(VirtualPosition),
    Done,
}

/// An iterator over records of a BAM reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'a, R>
where
    R: Read + Seek,
{
    reader: &'a mut Reader<bgzf::Reader<R>>,

    chunks: Vec<Chunk>,
    i: usize,

    reference_sequence_id: usize,
    interval: Interval,

    state: State,
    record: Record,
}

impl<'a, R> Query<'a, R>
where
    R: Read + Seek,
{
    pub(super) fn new(
        reader: &'a mut Reader<bgzf::Reader<R>>,
        chunks: Vec<Chunk>,
        reference_sequence_id: usize,
        interval: Interval,
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
                State::Read(chunk_end) => match self.next_record() {
                    Ok(Some(record)) => {
                        if self.reader.virtual_position() >= chunk_end {
                            self.state = State::Seek;
                        }

                        if intersects(&record, self.reference_sequence_id, self.interval) {
                            return Some(Ok(record));
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

pub(crate) fn intersects(
    record: &Record,
    reference_sequence_id: usize,
    interval: Interval,
) -> bool {
    match (
        record.reference_sequence_id(),
        record.alignment_start(),
        record.alignment_end(),
    ) {
        (Some(id), Some(start), Some(end)) => {
            id == reference_sequence_id && in_interval(start, end, interval)
        }
        _ => false,
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_in_interval() -> Result<(), noodles_core::position::TryFromIntError> {
        fn p(n: usize) -> Result<Position, noodles_core::position::TryFromIntError> {
            Position::try_from(n)
        }

        let start = p(5)?;
        let end = p(8)?;

        assert!(in_interval(start, end, p(4)?..p(7)?));
        assert!(in_interval(start, end, p(4)?..=p(6)?));

        assert!(in_interval(start, end, p(6)?..p(8)?));
        assert!(in_interval(start, end, p(6)?..=p(7)?));

        assert!(in_interval(start, end, p(7)?..p(10)?));
        assert!(in_interval(start, end, p(7)?..=p(9)?));

        assert!(in_interval(start, end, p(4)?..p(10)?));
        assert!(in_interval(start, end, p(4)?..=p(9)?));

        assert!(in_interval(start, end, p(4)?..));

        assert!(in_interval(start, end, ..p(6)?));
        assert!(in_interval(start, end, ..p(10)?));

        assert!(in_interval(start, end, ..=p(5)?));
        assert!(in_interval(start, end, ..=p(9)?));

        assert!(!in_interval(start, end, p(2)?..p(5)?));
        assert!(!in_interval(start, end, p(2)?..=p(4)?));

        assert!(!in_interval(start, end, p(9)?..p(12)?));
        assert!(!in_interval(start, end, p(9)?..=p(11)?));

        assert!(!in_interval(start, end, p(9)?..));

        assert!(!in_interval(start, end, ..p(5)?));
        assert!(!in_interval(start, end, ..=p(4)?));

        Ok(())
    }
}
