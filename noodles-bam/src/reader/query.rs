use std::io::{self, Read, Seek};

use noodles_bgzf::VirtualPosition;

use crate::{bai::index::reference::bin::Chunk, Record};

use super::Reader;

enum State {
    Seek,
    Read,
    End,
}

pub struct Query<'a, R: Read + Seek> {
    reader: &'a mut Reader<R>,
    chunks: Vec<Chunk>,
    reference_sequence_id: usize,
    start: u64,
    end: u64,
    i: usize,
    current_chunk: Chunk,
    state: State,
    record: Record,
}

impl<'a, R: Read + Seek> Query<'a, R> {
    pub fn new(
        reader: &'a mut Reader<R>,
        chunks: Vec<Chunk>,
        reference_sequence_id: usize,
        start: u64,
        end: u64,
    ) -> Self {
        let current_chunk = Chunk::new(VirtualPosition::from(0), VirtualPosition::from(1));

        Self {
            reader,
            chunks,
            reference_sequence_id,
            start,
            end,
            i: 0,
            current_chunk,
            state: State::Seek,
            record: Record::default(),
        }
    }

    fn next_chunk(&mut self) -> io::Result<bool> {
        if self.i >= self.chunks.len() {
            return Ok(false);
        }

        self.current_chunk = self.chunks[self.i];

        let pos = self.current_chunk.start();
        self.reader.seek(pos)?;

        self.i += 1;

        Ok(true)
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
                    self.state = match self.next_chunk() {
                        Ok(has_next_chunk) => {
                            if has_next_chunk {
                                State::Read
                            } else {
                                State::End
                            }
                        }
                        Err(e) => return Some(Err(e)),
                    }
                }
                State::Read => {
                    let result = self.read_record();

                    if self.reader.virtual_position() >= self.current_chunk.end() {
                        self.state = State::Seek;
                    }

                    if let Some(Ok(record)) = result {
                        let reference_sequence_id = record.reference_sequence_id() as usize;

                        let record_start = (record.position() + 1) as u64;
                        let record_mapped_len = record.cigar().mapped_len() as u64;
                        let record_end = record_start + record_mapped_len + 1;

                        if reference_sequence_id == self.reference_sequence_id
                            && in_interval(record_start, record_end, self.start, self.end)
                        {
                            return Some(Ok(record));
                        }
                    }
                }
                State::End => return None,
            }
        }
    }
}

fn in_interval(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> bool {
    a_start <= b_end && b_start <= a_end
}
