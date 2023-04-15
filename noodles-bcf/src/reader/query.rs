use std::{
    io::{self, Read, Seek},
    vec,
};

use noodles_bgzf as bgzf;
use noodles_core::{region::Interval, Position};
use noodles_csi::index::reference_sequence::bin::Chunk;
use noodles_vcf as vcf;

use crate::header::StringMaps;

use super::Reader;

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

/// An iterator over records of a BCF reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'r, 'h, R>
where
    R: Read + Seek,
{
    reader: &'r mut Reader<bgzf::Reader<R>>,

    header: &'h vcf::Header,
    chunks: vec::IntoIter<Chunk>,

    chromosome_id: usize,
    interval: Interval,

    state: State,
    record: vcf::Record,
}

impl<'r, 'h, R> Query<'r, 'h, R>
where
    R: Read + Seek,
{
    pub(crate) fn new(
        reader: &'r mut Reader<bgzf::Reader<R>>,
        header: &'h vcf::Header,
        chunks: Vec<Chunk>,
        chromosome_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader,

            header,
            chunks: chunks.into_iter(),

            chromosome_id,
            interval,

            state: State::Seek,
            record: vcf::Record::default(),
        }
    }

    fn read_record(&mut self) -> io::Result<Option<vcf::Record>> {
        self.reader
            .read_record(self.header, &mut self.record)
            .map(|n| match n {
                0 => None,
                _ => Some(self.record.clone()),
            })
    }
}

impl<'r, 'h, R> Iterator for Query<'r, 'h, R>
where
    R: Read + Seek,
{
    type Item = io::Result<vcf::Record>;

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

                        match intersects(
                            self.reader.string_maps(),
                            &record,
                            self.chromosome_id,
                            self.interval,
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

fn intersects(
    string_maps: &StringMaps,
    record: &vcf::Record,
    chromosome_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    let chromosome = record.chromosome();

    let id = string_maps
        .contigs()
        .get_index_of(&chromosome.to_string())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("chromosome does not exist in contigs: {chromosome}"),
            )
        })?;

    let start = Position::try_from(usize::from(record.position()))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let end = record
        .end()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .map(usize::from)
        .and_then(|n| {
            Position::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

    let record_interval = Interval::from(start..=end);

    Ok(id == chromosome_id && record_interval.intersects(region_interval))
}
