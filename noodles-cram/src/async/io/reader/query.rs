use std::{io::SeekFrom, slice, vec};

use futures::{Stream, stream};
use noodles_core::region::Interval;
use noodles_sam as sam;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::{crai, io::reader::Container};

/// An async reader over records of an async CRAM reader that intersects the given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'r, 'h: 'r, 'i: 'r, R> {
    inner: &'r mut Reader<R>,

    header: &'h sam::Header,

    index: slice::Iter<'i, crai::Record>,

    reference_sequence_id: usize,
    interval: Interval,

    records: vec::IntoIter<sam::alignment::RecordBuf>,
}

impl<'r, 'h: 'r, 'i: 'r, R> Query<'r, 'h, 'i, R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    pub(super) fn new(
        inner: &'r mut Reader<R>,
        header: &'h sam::Header,
        index: &'i crai::Index,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            inner,
            header,
            index: index.iter(),
            reference_sequence_id,
            interval,
            records: Vec::new().into_iter(),
        }
    }

    /// Reads a record.
    pub async fn read_record_buf(
        &mut self,
        record: &mut sam::alignment::RecordBuf,
    ) -> io::Result<usize> {
        loop {
            match self.records.next() {
                Some(r) => {
                    if let (Some(start), Some(end)) = (r.alignment_start(), r.alignment_end()) {
                        let alignment_interval = (start..=end).into();

                        if self.interval.intersects(alignment_interval) {
                            *record = r;
                            return Ok(1);
                        }
                    }
                }
                None => match read_next_container(self).await {
                    Some(Ok(())) => {}
                    Some(Err(e)) => return Err(e),
                    None => return Ok(0),
                },
            }
        }
    }

    /// Returns a stream over records.
    pub fn records(self) -> impl Stream<Item = io::Result<sam::alignment::RecordBuf>>
    where
        R: AsyncRead + AsyncSeek + Unpin,
    {
        Box::pin(stream::try_unfold(self, |mut ctx| async {
            let mut record = sam::alignment::RecordBuf::default();

            match ctx.read_record_buf(&mut record).await? {
                0 => Ok(None),
                _ => Ok(Some((record, ctx))),
            }
        }))
    }
}

async fn read_next_container<R>(ctx: &mut Query<'_, '_, '_, R>) -> Option<io::Result<()>>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    let index_record = ctx.index.next()?;

    if index_record.reference_sequence_id() != Some(ctx.reference_sequence_id) {
        return Some(Ok(()));
    }

    if let Err(e) = ctx.inner.seek(SeekFrom::Start(index_record.offset())).await {
        return Some(Err(e));
    }

    let mut container = Container::default();

    match ctx.inner.read_container(&mut container).await {
        Ok(0) => return None,
        Ok(_) => {}
        Err(e) => return Some(Err(e)),
    };

    let compression_header = match container.compression_header() {
        Ok(compression_header) => compression_header,
        Err(e) => return Some(Err(e)),
    };

    let records = container
        .slices()
        .map(|result| {
            let slice = result?;

            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

            slice
                .records(
                    ctx.inner.reference_sequence_repository.clone(),
                    ctx.header,
                    &compression_header,
                    &core_data_src,
                    &external_data_srcs,
                )
                .and_then(|records| {
                    records
                        .into_iter()
                        .map(|record| {
                            sam::alignment::RecordBuf::try_from_alignment_record(
                                ctx.header, &record,
                            )
                        })
                        .collect::<io::Result<Vec<_>>>()
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
