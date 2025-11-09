use futures::{Stream, stream};
use noodles_core::Region;
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, Lines};

use crate::{
    binning_index::index::{Header, header::format::CoordinateSystem},
    io::{
        IndexedRecord,
        indexed_records::{Record, parse_record},
    },
};

struct Context<'r, R> {
    lines: Lines<R>,
    line_comment_prefix: char,
    reference_sequence_name_index: usize,
    start_position_index: usize,
    end_position_index: Option<usize>,
    coordinate_system: CoordinateSystem,
    region: &'r Region,
}

pub(super) fn filtered_indexed_records<R>(
    reader: R,
    header: &Header,
    region: &Region,
) -> impl Stream<Item = io::Result<Record>>
where
    R: AsyncBufRead + Unpin,
{
    let ctx = Context {
        lines: reader.lines(),
        line_comment_prefix: char::from(header.line_comment_prefix()),
        reference_sequence_name_index: header.reference_sequence_name_index(),
        start_position_index: header.start_position_index(),
        end_position_index: header.end_position_index(),
        coordinate_system: header.format().coordinate_system(),
        region,
    };

    Box::pin(stream::try_unfold(ctx, |mut ctx| async move {
        loop {
            let line = match ctx.lines.next_line().await? {
                Some(s) => s,
                None => return Ok(None),
            };

            if line.starts_with(ctx.line_comment_prefix) {
                continue;
            }

            let record = match parse_record(
                line,
                ctx.reference_sequence_name_index,
                ctx.start_position_index,
                ctx.end_position_index,
                ctx.coordinate_system,
            ) {
                Ok(record) => record,
                Err(e) => return Err(io::Error::new(io::ErrorKind::InvalidData, e)),
            };

            if !intersects(&record, ctx.region) {
                continue;
            }

            return Ok(Some((record, ctx)));
        }
    }))
}

fn intersects(record: &Record, region: &Region) -> bool {
    record.indexed_reference_sequence_name() == region.name()
        && record.indexed_interval().intersects(region.interval())
}
