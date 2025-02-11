use std::{io::SeekFrom, slice, vec};

use futures::{stream, Stream};
use noodles_core::region::Interval;
use noodles_sam as sam;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::{crai, Record};

struct Context<'r, 'h: 'r, 'i: 'r, R> {
    reader: &'r mut Reader<R>,

    header: &'h sam::Header,

    index: slice::Iter<'i, crai::Record>,

    reference_sequence_id: usize,
    interval: Interval,

    records: vec::IntoIter<Record>,
}

pub(super) fn query<'r, 'h: 'r, 'i: 'r, R>(
    reader: &'r mut Reader<R>,
    header: &'h sam::Header,
    index: &'i crai::Index,
    reference_sequence_id: usize,
    interval: Interval,
) -> impl Stream<Item = io::Result<Record>> + 'r
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
                    if let (Some(start), Some(end)) = (r.alignment_start(), r.alignment_end()) {
                        let alignment_interval = (start..=end).into();

                        if ctx.interval.intersects(alignment_interval) {
                            return Ok(Some((r, ctx)));
                        }
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

    if index_record.reference_sequence_id() != Some(ctx.reference_sequence_id) {
        return Some(Ok(()));
    }

    if let Err(e) = ctx
        .reader
        .seek(SeekFrom::Start(index_record.offset()))
        .await
    {
        return Some(Err(e));
    }

    let container = match ctx.reader.read_container().await {
        Ok(Some(c)) => c,
        Ok(None) => return None,
        Err(e) => return Some(Err(e)),
    };

    let records = container
        .slices()
        .iter()
        .map(|slice| {
            let compression_header = container.compression_header();

            slice.records(compression_header).and_then(|mut records| {
                slice.resolve_records(
                    ctx.reader.reference_sequence_repository(),
                    ctx.header,
                    compression_header,
                    &mut records,
                )?;

                Ok(records)
            })
        })
        .collect::<Result<Vec<_>, _>>();

    let records = match records {
        Ok(records) => records,
        Err(e) => return Some(Err(e)),
    };

    ctx.records = records
        .into_iter()
        .flatten()
        .collect::<Vec<_>>()
        .into_iter();

    Some(Ok(()))
}
