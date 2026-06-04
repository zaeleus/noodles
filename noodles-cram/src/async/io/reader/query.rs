use std::{io::SeekFrom, slice, vec};

use futures::{Stream, stream};
use noodles_core::region::Interval;
use noodles_sam as sam;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::{crai, io::reader::Container};

struct Context<'r, 'h: 'r, 'i: 'r, R> {
    reader: &'r mut Reader<R>,

    header: &'h sam::Header,

    index: slice::Iter<'i, crai::Record>,

    reference_sequence_id: usize,
    interval: Interval,

    records: vec::IntoIter<sam::alignment::RecordBuf>,
}

pub(super) fn query<'r, 'h: 'r, 'i: 'r, R>(
    reader: &'r mut Reader<R>,
    header: &'h sam::Header,
    index: &'i crai::Index,
    reference_sequence_id: usize,
    interval: Interval,
) -> impl Stream<Item = io::Result<sam::alignment::RecordBuf>> + 'r
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    let ctx = Context {
        reader,

        header,

        index: index.iter(),

        reference_sequence_id,
        interval,

        records: Vec::new().into_iter(),
    };

    Box::pin(stream::try_unfold(ctx, |mut ctx| async {
        loop {
            match ctx.records.next() {
                Some(r) => {
                    if intersects(&r, ctx.reference_sequence_id, ctx.interval) {
                        return Ok(Some((r, ctx)));
                    }
                }
                None => match read_next_container(&mut ctx).await {
                    Some(Ok(())) => {}
                    Some(Err(e)) => return Err(e),
                    None => return Ok(None),
                },
            }
        }
    }))
}

async fn read_next_container<R>(ctx: &mut Context<'_, '_, '_, R>) -> Option<io::Result<()>>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    let index_record = ctx.index.next()?;

    if !intersects_index_record(index_record, ctx.reference_sequence_id, ctx.interval) {
        return Some(Ok(()));
    }

    if let Err(e) = ctx
        .reader
        .seek(SeekFrom::Start(index_record.offset()))
        .await
    {
        return Some(Err(e));
    }

    let mut container = Container::default();

    match ctx.reader.read_container(&mut container).await {
        Ok(0) => return None,
        Ok(_) => {}
        Err(e) => return Some(Err(e)),
    };

    let compression_header = match container.compression_header() {
        Ok(compression_header) => compression_header,
        Err(e) => return Some(Err(e)),
    };

    let slice = match container.slice(index_record.landmark()) {
        Ok(slice) => slice,
        Err(e) => return Some(Err(e)),
    };

    let (core_data_src, external_data_srcs) = match slice.decode_blocks() {
        Ok(blocks) => blocks,
        Err(e) => return Some(Err(e)),
    };

    let records = slice
        .records(
            ctx.reader.reference_sequence_repository.clone(),
            ctx.header,
            &compression_header,
            &core_data_src,
            &external_data_srcs,
        )
        .and_then(|records| {
            records
                .into_iter()
                .map(|record| {
                    sam::alignment::RecordBuf::try_from_alignment_record(ctx.header, &record)
                })
                .collect::<io::Result<Vec<_>>>()
        });

    let records = match records {
        Ok(records) => records,
        Err(e) => return Some(Err(e)),
    };

    ctx.records = records.into_iter();

    Some(Ok(()))
}

fn intersects(
    record: &sam::alignment::RecordBuf,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> bool {
    match (
        record.reference_sequence_id(),
        record.alignment_start(),
        record.alignment_end(),
    ) {
        (Some(id), Some(start), Some(end)) if id == reference_sequence_id => {
            let alignment_interval = (start..=end).into();
            region_interval.intersects(alignment_interval)
        }
        _ => false,
    }
}

fn intersects_index_record(
    record: &crai::Record,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> bool {
    if record.reference_sequence_id() != Some(reference_sequence_id) {
        return false;
    }

    let Some(start) = record.alignment_start() else {
        return false;
    };

    let Some(end) = start.checked_add(record.alignment_span().saturating_sub(1)) else {
        return true;
    };

    region_interval.intersects((start..=end).into())
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;
    use noodles_sam::alignment::{
        record::cigar::{Op, op::Kind},
        record_buf::Cigar,
    };

    use super::*;

    #[test]
    fn test_intersects() {
        let record = build_record(0, 8, 5);
        let interval =
            Interval::from(Position::try_from(10).unwrap()..=Position::try_from(13).unwrap());

        assert!(intersects(&record, 0, interval));
    }

    #[test]
    fn test_intersects_with_different_reference_sequence() {
        let record = build_record(1, 8, 5);
        let interval =
            Interval::from(Position::try_from(10).unwrap()..=Position::try_from(13).unwrap());

        assert!(!intersects(&record, 0, interval));
    }

    #[test]
    fn test_intersects_index_record() {
        let record = crai::Record::new(Some(0), Position::new(8), 5, 13, 21, 34);
        let interval =
            Interval::from(Position::try_from(10).unwrap()..=Position::try_from(13).unwrap());

        assert!(intersects_index_record(&record, 0, interval));
    }

    #[test]
    fn test_intersects_index_record_with_nonoverlapping_interval() {
        let record = crai::Record::new(Some(0), Position::new(8), 5, 13, 21, 34);
        let interval =
            Interval::from(Position::try_from(13).unwrap()..=Position::try_from(21).unwrap());

        assert!(!intersects_index_record(&record, 0, interval));
    }

    #[test]
    fn test_intersects_index_record_with_different_reference_sequence() {
        let record = crai::Record::new(Some(1), Position::new(8), 5, 13, 21, 34);
        let interval = Interval::from(Position::MIN..=Position::MAX);

        assert!(!intersects_index_record(&record, 0, interval));
    }

    fn build_record(
        reference_sequence_id: usize,
        alignment_start: usize,
        alignment_span: usize,
    ) -> sam::alignment::RecordBuf {
        let cigar: Cigar = [Op::new(Kind::Match, alignment_span)].into_iter().collect();

        sam::alignment::RecordBuf::builder()
            .set_reference_sequence_id(reference_sequence_id)
            .set_alignment_start(Position::try_from(alignment_start).unwrap())
            .set_cigar(cigar)
            .build()
    }
}
