use std::{
    io::{self, Read, Seek},
    ops::{Bound, RangeBounds},
};

use noodles_bgzf as bgzf;
use noodles_csi::index::reference_sequence::bin::Chunk;

use super::Reader;
use crate::{Header, Record};

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

/// An iterator over records of a VCF reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'r, 'h, R>
where
    R: Read + Seek + 'r,
{
    reader: &'r mut Reader<bgzf::Reader<R>>,

    chunks: Vec<Chunk>,
    i: usize,

    reference_sequence_name: String,
    start: i32,
    end: i32,

    state: State,
    header: &'h Header,
    line_buf: String,
}

impl<'r, 'h, R> Query<'r, 'h, R>
where
    R: Read + Seek,
{
    pub(super) fn new<B>(
        reader: &'r mut Reader<bgzf::Reader<R>>,
        chunks: Vec<Chunk>,
        reference_sequence_name: String,
        interval: B,
        header: &'h Header,
    ) -> Self
    where
        B: RangeBounds<i32>,
    {
        let (start, end) = resolve_interval(interval);

        Self {
            reader,

            chunks,
            i: 0,

            reference_sequence_name,
            start,
            end,

            state: State::Seek,
            header,
            line_buf: String::new(),
        }
    }

    fn read_record(&mut self) -> io::Result<Option<Record>> {
        self.line_buf.clear();

        self.reader
            .read_record(&mut self.line_buf)
            .and_then(|n| match n {
                0 => Ok(None),
                _ => Record::try_from_str(&self.line_buf, self.header)
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }
}

impl<'r, 'h, R> Iterator for Query<'r, 'h, R>
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

                        match intersects(
                            &record,
                            &self.reference_sequence_name,
                            self.start,
                            self.end,
                        ) {
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

pub(crate) fn resolve_interval<B>(interval: B) -> (i32, i32)
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

pub(crate) fn next_chunk(chunks: &[Chunk], i: &mut usize) -> Option<Chunk> {
    let chunk = chunks.get(*i).copied();
    *i += 1;
    chunk
}

pub(crate) fn intersects(
    record: &Record,
    reference_sequence_name: &str,
    interval_start: i32,
    interval_end: i32,
) -> io::Result<bool> {
    let name = record.chromosome().to_string();

    let start = i32::from(record.position());
    let end = record
        .end()
        .map(i32::from)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(name == reference_sequence_name && in_interval(start, end, interval_start, interval_end))
}

fn in_interval(a_start: i32, a_end: i32, b_start: i32, b_end: i32) -> bool {
    a_start <= b_end && b_start <= a_end
}
