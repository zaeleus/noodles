use std::io::{self, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_tabix::index::reference_sequence::bin::Chunk;

use crate::{record::Chromosome, Record};

use super::Reader;

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    End,
}

/// An iterator over records of a VCF reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'a, R>
where
    R: Read + Seek,
{
    reader: &'a mut Reader<bgzf::Reader<R>>,
    chunks: Vec<Chunk>,
    reference_sequence_name: String,
    start: i32,
    end: i32,
    i: usize,
    state: State,
    line_buf: String,
}

impl<'a, R> Query<'a, R>
where
    R: Read + Seek,
{
    pub(crate) fn new(
        reader: &'a mut Reader<bgzf::Reader<R>>,
        chunks: Vec<Chunk>,
        reference_sequence_name: String,
        start: i32,
        end: i32,
    ) -> Self {
        Self {
            reader,
            chunks,
            reference_sequence_name,
            start,
            end,
            i: 0,
            state: State::Seek,
            line_buf: String::new(),
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

    fn read_and_parse_record(&mut self) -> Option<io::Result<Record>> {
        self.line_buf.clear();

        match self.reader.read_record(&mut self.line_buf) {
            Ok(0) => None,
            Ok(_) => Some(
                self.line_buf
                    .parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            ),
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
                State::Read(chunk_end) => match self.read_and_parse_record() {
                    Some(result) => {
                        if self.reader.virtual_position() >= chunk_end {
                            self.state = State::Seek;
                        }

                        match result {
                            Ok(record) => {
                                let reference_sequence_name = match record.chromosome() {
                                    Chromosome::Name(n) => n.into(),
                                    Chromosome::Symbol(n) => n.to_string(),
                                };

                                let start = i32::from(record.position());
                                let end = start + 1; // TODO: records with END

                                if reference_sequence_name == self.reference_sequence_name
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
