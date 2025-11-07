//! Queries an indexed TSV with the given region.
//!
//! The input must have an associated tabix index in the same directory.
//!
//! The result matches the output of `tabix <src> <region>`.

use std::env;

use noodles_core::{Position, Region, region::Interval};
use noodles_csi::{
    self as csi, BinningIndex, binning_index::index::header::format::CoordinateSystem,
};
use noodles_tabix as tabix;
use tokio::{
    fs::File,
    io::{self, AsyncBufReadExt},
};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");
    let region = args.next().expect("missing region").parse()?;

    let index_src = format!("{src}.tbi");
    let index = tabix::fs::read(index_src)?;
    let header = index.header().expect("missing index header").clone();

    let mut reader = File::open(src)
        .await
        .map(|file| csi::r#async::io::IndexedReader::new(file, index))?;

    let line_comment_prefix = char::from(header.line_comment_prefix());

    let query = reader.query(&region)?;
    let mut lines = query.lines();

    while let Some(line) = lines.next_line().await? {
        if line.starts_with(line_comment_prefix) {
            continue;
        }

        let (reference_sequence_name, interval) = parse_record(&line, &header)?;

        if intersects(reference_sequence_name, interval, &region) {
            println!("{line}");
        }
    }

    Ok(())
}

fn parse_record<'r>(
    line: &'r str,
    header: &csi::binning_index::index::Header,
) -> io::Result<(&'r str, Interval)> {
    const SEPARATOR: char = '\t';

    let fields: Vec<_> = line.split(SEPARATOR).collect();
    let reference_sequence_name = &fields[header.reference_sequence_name_index()];

    let mut start = parse_position(fields[header.start_position_index()])?;

    if header.format().coordinate_system() == CoordinateSystem::Bed {
        start = start.checked_add(1).expect("attempt to add with overflow")
    }

    let end = if let Some(i) = header.end_position_index() {
        parse_position(fields[i])?
    } else {
        start
    };

    Ok((reference_sequence_name, (start..=end).into()))
}

fn parse_position(s: &str) -> io::Result<Position> {
    s.parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
}

fn intersects(reference_sequence_name: &str, interval: Interval, region: &Region) -> bool {
    reference_sequence_name == region.name() && interval.intersects(region.interval())
}
