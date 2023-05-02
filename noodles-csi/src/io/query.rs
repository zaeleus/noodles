use std::{
    io::{self, BufRead, Read, Seek},
    vec,
};

use noodles_bgzf as bgzf;

use crate::index::reference_sequence::bin::Chunk;

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

/// A query reader.
///
/// This reader returns the uncompressed data between all the given chunks.
pub struct Query<'r, R> {
    reader: &'r mut bgzf::Reader<R>,
    chunks: vec::IntoIter<Chunk>,
    state: State,
}

impl<'r, R> Query<'r, R>
where
    R: Read + Seek,
{
    /// Creates a query reader.
    pub fn new(reader: &'r mut bgzf::Reader<R>, chunks: Vec<Chunk>) -> Self {
        Self {
            reader,
            chunks: chunks.into_iter(),
            state: State::Seek,
        }
    }
}

impl<'r, R> Read for Query<'r, R>
where
    R: Read + Seek,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut src = self.fill_buf()?;
        let amt = src.read(buf)?;
        self.consume(amt);
        Ok(amt)
    }
}

impl<'r, R> BufRead for Query<'r, R>
where
    R: Read + Seek,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        loop {
            match self.state {
                State::Seek => {
                    self.state = match self.chunks.next() {
                        Some(chunk) => {
                            self.reader.seek(chunk.start())?;
                            State::Read(chunk.end())
                        }
                        None => State::Done,
                    }
                }
                State::Read(chunk_end) => {
                    if self.reader.virtual_position() < chunk_end {
                        return self.reader.fill_buf();
                    } else {
                        self.state = State::Seek;
                    }
                }
                State::Done => return Ok(&[]),
            }
        }
    }

    fn consume(&mut self, amt: usize) {
        self.reader.consume(amt);
    }
}
