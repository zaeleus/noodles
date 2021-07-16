use std::{
    io::{self, Read, Seek},
    ops::{Bound, RangeBounds},
};

use noodles_bgzf::VirtualPosition;
use noodles_csi::index::reference_sequence::bin::Chunk;

use crate::Record;

use super::Reader;

enum State {
    Seek,
    Read(VirtualPosition),
    End,
}

/// An iterator over records of a BAM reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'a, R>
where
    R: Read + Seek,
{
    reader: &'a mut Reader<R>,
    chunks: Vec<Chunk>,
    reference_sequence_id: usize,
    start: i32,
    end: i32,
    i: usize,
    state: State,
    record: Record,
}

impl<'a, R> Query<'a, R>
where
    R: Read + Seek,
{
    pub(crate) fn new<B>(
        reader: &'a mut Reader<R>,
        chunks: Vec<Chunk>,
        reference_sequence_id: usize,
        interval: B,
    ) -> Self
    where
        B: RangeBounds<i32>,
    {
        let (start, end) = match (interval.start_bound(), interval.end_bound()) {
            (Bound::Included(s), Bound::Included(e)) => (*s, *e),
            (Bound::Included(s), Bound::Unbounded) => (*s, i32::MAX),
            (Bound::Unbounded, Bound::Unbounded) => (1, i32::MAX),
            _ => todo!(),
        };

        Self {
            reader,
            chunks,
            reference_sequence_id,
            start,
            end,
            i: 0,
            state: State::Seek,
            record: Record::default(),
        }
    }

    fn next_chunk(&mut self) -> io::Result<Option<VirtualPosition>> {
        if self.i >= self.chunks.len() {
            return Ok(None);
        }

        let chunk = self.chunks[self.i];
        self.reader.seek(chunk.start())?;

        self.i += 1;

        Ok(Some(chunk.end()))
    }

    fn read_record(&mut self) -> Option<io::Result<Record>> {
        match self.reader.read_record(&mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
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
                    self.state = match self.next_chunk() {
                        Ok(Some(chunk_end)) => State::Read(chunk_end),
                        Ok(None) => State::End,
                        Err(e) => return Some(Err(e)),
                    }
                }
                State::Read(chunk_end) => match self.read_record() {
                    Some(result) => {
                        if self.reader.virtual_position() >= chunk_end {
                            self.state = State::Seek;
                        }

                        match result {
                            Ok(record) => {
                                if let Some(reference_sequence_id) = record.reference_sequence_id()
                                {
                                    let reference_sequence_id =
                                        i32::from(reference_sequence_id) as usize;

                                    let record_start =
                                        record.position().map(i32::from).expect("missing position");
                                    let record_reference_len = match record.cigar().reference_len()
                                    {
                                        Ok(len) => len as i32,
                                        Err(e) => return Some(Err(e)),
                                    };
                                    let record_end = record_start + record_reference_len - 1;

                                    if reference_sequence_id == self.reference_sequence_id
                                        && in_interval(
                                            record_start,
                                            record_end,
                                            self.start,
                                            self.end,
                                        )
                                    {
                                        return Some(Ok(record));
                                    }
                                }
                            }
                            Err(e) => return Some(Err(e)),
                        }
                    }
                    None => {
                        self.state = State::Seek;
                    }
                },
                State::End => return None,
            }
        }
    }
}

fn in_interval(a_start: i32, a_end: i32, b_start: i32, b_end: i32) -> bool {
    a_start <= b_end && b_start <= a_end
}
