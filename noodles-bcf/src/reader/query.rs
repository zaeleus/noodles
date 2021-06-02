use std::io::{self, Read, Seek};

use noodles_bgzf::{self as bgzf, index::Chunk};

use crate::Record;

use super::Reader;

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    End,
}

/// An iterator over records of a BCF reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'a, R>
where
    R: Read + Seek,
{
    reader: &'a mut Reader<R>,
    chunks: Vec<Chunk>,
    chromosome_id: usize,
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
    pub(crate) fn new(
        reader: &'a mut Reader<R>,
        chunks: Vec<Chunk>,
        chromosome_id: usize,
        start: i32,
        end: i32,
    ) -> Self {
        Self {
            reader,
            chunks,
            chromosome_id,
            start,
            end,
            i: 0,
            state: State::Seek,
            record: Record::default(),
        }
    }

    fn next_chunk(&mut self) -> io::Result<Option<bgzf::VirtualPosition>> {
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
                                let chromosome_id = match record.chromosome_id() {
                                    Ok(id) => id as usize,
                                    Err(e) => return Some(Err(e)),
                                };

                                let start = match record.position() {
                                    Ok(pos) => i32::from(pos),
                                    Err(e) => return Some(Err(e)),
                                };

                                let end = match record.end() {
                                    Ok(pos) => i32::from(pos),
                                    Err(e) => return Some(Err(e)),
                                };

                                if chromosome_id == self.chromosome_id
                                    && in_interval(start, end, self.start, self.end)
                                {
                                    return Some(Ok(record));
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
