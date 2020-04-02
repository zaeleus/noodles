use std::io::{self, Read, Seek};

pub use crate::Record;

use crate::bai;

use super::Reader;

enum State {
    Seek,
    Read,
    End,
}

pub struct Query<'a, R: Read + Seek> {
    reader: &'a mut Reader<R>,
    chunks: Vec<bai::Chunk>,
    i: usize,
    current_chunk: bai::Chunk,
    state: State,
    record: Record,
}

impl<'a, R: Read + Seek> Query<'a, R> {
    pub fn new(reader: &'a mut Reader<R>, chunks: Vec<bai::Chunk>) -> Self {
        let current_chunk = bai::Chunk::new(0, 1);

        Self {
            reader,
            chunks,
            i: 0,
            current_chunk,
            state: State::Seek,
            record: Record::new(),
        }
    }

    fn next_chunk(&mut self) -> bool {
        if self.i >= self.chunks.len() {
            return false;
        }

        self.current_chunk = self.chunks[self.i];

        let pos = self.current_chunk.start();
        self.reader.seek(pos).unwrap();

        self.i += 1;

        true
    }

    fn read_record(&mut self) -> Option<io::Result<Record>> {
        match self.reader.read_record(&mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}

impl<'a, R: Read + Seek> Iterator for Query<'a, R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.state {
                State::Seek => {
                    self.state = if self.next_chunk() {
                        State::Read
                    } else {
                        State::End
                    };
                }
                State::Read => {
                    let result = self.read_record();

                    if self.reader.virtual_position() >= self.current_chunk.end() {
                        self.state = State::Seek;
                    }

                    return result;
                }
                State::End => return None,
            }
        }
    }
}
