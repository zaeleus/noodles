use std::{
    io::{self, Read, Seek},
    vec,
};

use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
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

    chunks: vec::IntoIter<Chunk>,

    reference_sequence_name: String,
    interval: Interval,

    state: State,
    header: &'h Header,
    line_buf: String,
}

impl<'r, 'h, R> Query<'r, 'h, R>
where
    R: Read + Seek,
{
    pub(super) fn new(
        reader: &'r mut Reader<bgzf::Reader<R>>,
        chunks: Vec<Chunk>,
        reference_sequence_name: String,
        interval: Interval,
        header: &'h Header,
    ) -> Self {
        Self {
            reader,

            chunks: chunks.into_iter(),

            reference_sequence_name,
            interval,

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
                    self.state = match self.chunks.next() {
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

                        match intersects(&record, &self.reference_sequence_name, self.interval) {
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

pub(crate) fn intersects(
    record: &Record,
    reference_sequence_name: &str,
    region_interval: Interval,
) -> io::Result<bool> {
    use noodles_core::Position;

    let name = record.chromosome().to_string();

    let start = Position::try_from(usize::from(record.position()))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let end = record
        .end()
        .map(usize::from)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|n| {
            Position::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

    let record_interval = Interval::from(start..=end);

    Ok(name == reference_sequence_name && record_interval.intersects(region_interval))
}
