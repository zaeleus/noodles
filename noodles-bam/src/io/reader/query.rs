use std::io;

use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_sam::alignment::Record as _;

use crate::Record;

struct Reader<'r, R> {
    inner: super::Reader<csi::io::Query<'r, R>>,
    reference_sequence_id: usize,
    interval: Interval,
}

impl<'r, R> Reader<'r, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    fn new(
        reader: &'r mut R,
        chunks: Vec<Chunk>,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            inner: super::Reader::from(csi::io::Query::new(reader, chunks)),
            reference_sequence_id,
            interval,
        }
    }

    fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        next_record(
            &mut self.inner,
            record,
            self.reference_sequence_id,
            self.interval,
        )
    }
}

/// An iterator over records of a BAM reader that intersects a given region.
///
/// This is created by calling [`super::Reader::query`].
pub struct Query<'r, R> {
    reader: Reader<'r, R>,
    record: Record,
}

impl<'r, R> Query<'r, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    pub(super) fn new(
        reader: &'r mut R,
        chunks: Vec<Chunk>,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader: Reader::new(reader, chunks, reference_sequence_id, interval),
            record: Record::default(),
        }
    }
}

impl<R> Iterator for Query<'_, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record(&mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}

pub(crate) fn intersects(
    record: &Record,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    let Some(id) = record.reference_sequence_id().transpose()? else {
        return Ok(false);
    };

    if id != reference_sequence_id {
        return Ok(false);
    }

    if interval_is_unbounded(region_interval) {
        Ok(true)
    } else {
        match (
            record.alignment_start().transpose()?,
            record.alignment_end().transpose()?,
        ) {
            (Some(start), Some(end)) => {
                let alignment_interval = (start..=end).into();
                Ok(region_interval.intersects(alignment_interval))
            }
            _ => Ok(false),
        }
    }
}

fn interval_is_unbounded(interval: Interval) -> bool {
    interval.start().is_none() && interval.end().is_none()
}

fn next_record<R>(
    reader: &mut super::Reader<csi::io::Query<'_, R>>,
    record: &mut Record,
    reference_sequence_id: usize,
    interval: Interval,
) -> io::Result<usize>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    loop {
        match reader.read_record(record)? {
            0 => return Ok(0),
            n => {
                if intersects(record, reference_sequence_id, interval)? {
                    return Ok(n);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{io::Cursor, num::NonZero};

    use noodles_core::Position;
    use noodles_csi::binning_index::Indexer;
    use noodles_sam::{
        self as sam,
        alignment::{
            RecordBuf,
            io::Write,
            record::{
                Flags,
                cigar::{Op, op::Kind},
            },
        },
        header::record::value::{Map, map::ReferenceSequence},
    };

    use super::*;
    use crate::{bai, io::Writer};

    fn write(header: &sam::Header, records: &[RecordBuf]) -> io::Result<Vec<u8>> {
        let mut writer = Writer::new(Vec::new());
        writer.write_header(header)?;

        for record in records {
            writer.write_alignment_record(header, record)?;
        }

        writer.into_inner().finish()
    }

    fn index(src: &[u8]) -> io::Result<bai::Index> {
        let mut reader = crate::io::Reader::new(src);
        let header = reader.read_header()?;

        let mut indexer = Indexer::default();
        let mut chunk_start = reader.get_ref().virtual_position();

        let mut record = Record::default();

        while reader.read_record(&mut record)? != 0 {
            let chunk_end = reader.get_ref().virtual_position();

            let alignment_context = match (
                record.reference_sequence_id().transpose()?,
                record.alignment_start().transpose()?,
                record.alignment_end().transpose()?,
            ) {
                (Some(id), Some(start), Some(end)) => {
                    let is_mapped = !record.flags().is_unmapped();
                    Some((id, start, end, is_mapped))
                }
                _ => None,
            };

            let chunk = Chunk::new(chunk_start, chunk_end);
            indexer.add_record(alignment_context, chunk)?;

            chunk_start = chunk_end;
        }

        Ok(indexer.build(header.reference_sequences().len()))
    }

    #[test]
    fn test_next() -> Result<(), Box<dyn std::error::Error>> {
        let header = sam::Header::builder()
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
            )
            .add_reference_sequence(
                "sq1",
                Map::<ReferenceSequence>::new(const { NonZero::new(13).unwrap() }),
            )
            .build();

        let records = [
            RecordBuf::builder()
                .set_reference_sequence_id(0)
                .set_flags(Flags::default())
                .set_alignment_start(Position::MIN)
                .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
                .build(),
            RecordBuf::builder()
                .set_reference_sequence_id(1)
                .set_flags(Flags::default())
                .set_alignment_start(Position::MIN)
                .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
                .build(),
            RecordBuf::builder()
                .set_reference_sequence_id(1)
                .set_flags(Flags::default())
                .set_alignment_start(Position::try_from(8)?)
                .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
                .build(),
        ];

        let src = write(&header, &records)?;
        let index = index(&src)?;

        let mut reader = crate::io::Reader::new(Cursor::new(src));

        let region = "sq1:2-5".parse()?;
        let query = reader.query(&header, &index, &region)?;

        let actual: Vec<_> = query
            .map(|result| {
                result.and_then(|record| RecordBuf::try_from_alignment_record(&header, &record))
            })
            .collect::<Result<_, _>>()?;

        let expected = [records[1].clone()];

        assert_eq!(actual, expected);

        Ok(())
    }
}
