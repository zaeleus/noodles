use std::vec;

use futures::{stream, Stream};
use noodles_sam as sam;
use tokio::io::{self, AsyncRead};

use super::Reader;
use crate::io::reader::Container;

struct Context<'r, 'h: 'r, R>
where
    R: AsyncRead + Unpin,
{
    reader: &'r mut Reader<R>,
    header: &'h sam::Header,
    container: Container,
    records: vec::IntoIter<sam::alignment::RecordBuf>,
}

pub fn records<'r, 'h: 'r, R>(
    reader: &'r mut Reader<R>,
    header: &'h sam::Header,
) -> impl Stream<Item = io::Result<sam::alignment::RecordBuf>> + 'r
where
    R: AsyncRead + Unpin,
{
    let ctx = Context {
        reader,
        header,
        container: Container::default(),
        records: Vec::new().into_iter(),
    };

    Box::pin(stream::try_unfold(ctx, |mut ctx| async {
        loop {
            match ctx.records.next() {
                Some(record) => return Ok(Some((record, ctx))),
                None => match read_next_container(&mut ctx).await {
                    Some(Ok(records)) => ctx.records = records.into_iter(),
                    Some(Err(e)) => return Err(e),
                    None => return Ok(None),
                },
            }
        }
    }))
}

async fn read_next_container<R>(
    ctx: &mut Context<'_, '_, R>,
) -> Option<io::Result<Vec<sam::alignment::RecordBuf>>>
where
    R: AsyncRead + Unpin,
{
    match ctx.reader.read_container(&mut ctx.container).await {
        Ok(0) => return None,
        Ok(_) => {}
        Err(e) => return Some(Err(e)),
    };

    let compression_header = match ctx.container.compression_header() {
        Ok(compression_header) => compression_header,
        Err(e) => return Some(Err(e)),
    };

    let records = ctx
        .container
        .slices()
        .map(|result| {
            let slice = result?;

            slice.records(&compression_header).and_then(|mut records| {
                slice.resolve_records(
                    ctx.reader.reference_sequence_repository(),
                    ctx.header,
                    &compression_header,
                    &mut records,
                )?;

                records
                    .into_iter()
                    .map(|record| {
                        sam::alignment::RecordBuf::try_from_alignment_record(ctx.header, &record)
                    })
                    .collect::<io::Result<Vec<_>>>()
            })
        })
        .collect::<Result<Vec<_>, _>>();

    let records = match records {
        Ok(records) => records.into_iter().flatten().collect::<Vec<_>>(),
        Err(e) => return Some(Err(e)),
    };

    Some(Ok(records))
}
