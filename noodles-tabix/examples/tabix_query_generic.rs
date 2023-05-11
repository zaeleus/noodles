//! Queries an indexed TSV with the given region.
//!
//! The input must have an associated tabix index in the same directory.
//!
//! The result matches the output of `tabix <src> <region>`.

use std::{
    env,
    fs::File,
    io::{self, BufRead},
};

use noodles_bgzf as bgzf;
use noodles_core::{region::Interval, Position, Region};
use noodles_csi::{self as csi, index::header::format::CoordinateSystem};
use noodles_tabix as tabix;

fn resolve_region(header: &csi::index::Header, region: &Region) -> io::Result<usize> {
    header
        .reference_sequence_names()
        .get_index_of(region.name())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "missing reference sequence name",
            )
        })
}

fn parse_start_position(s: &str, coordinate_system: CoordinateSystem) -> io::Result<Position> {
    fn invalid_position<E>(_: E) -> io::Error {
        io::Error::new(io::ErrorKind::InvalidData, "invalid position")
    }

    match coordinate_system {
        CoordinateSystem::Gff => s.parse().map_err(invalid_position),
        CoordinateSystem::Bed => s
            .parse::<usize>()
            .map_err(invalid_position)
            .and_then(|n| Position::try_from(n + 1).map_err(invalid_position)),
    }
}

fn intersects(
    header: &csi::index::Header,
    line: &str,
    region: &Region,
) -> Result<bool, Box<dyn std::error::Error>> {
    const DELIMITER: char = '\t';

    let fields: Vec<_> = line.split(DELIMITER).collect();

    let reference_sequence_name = fields[header.reference_sequence_name_index() - 1];

    let raw_start = fields[header.start_position_index() - 1];
    let coordinate_system = header.format().coordinate_system();
    let start = parse_start_position(raw_start, coordinate_system)?;

    let end = if let Some(i) = header.end_position_index() {
        fields[i - 1].parse()?
    } else {
        start.checked_add(1).expect("attempt to add with overflow")
    };

    let interval = Interval::from(start..=end);

    Ok(reference_sequence_name == region.name() && interval.intersects(region.interval()))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");
    let region = args.next().expect("missing region").parse()?;

    let tabix_src = format!("{src}.tbi");
    let index = tabix::read(tabix_src)?;
    let header = index
        .header()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing tabix header"))?;

    let mut reader = File::open(src).map(bgzf::Reader::new)?;

    let reference_sequence_id = resolve_region(header, &region)?;
    let chunks = index.query(reference_sequence_id, region.interval())?;
    let query = csi::io::Query::new(&mut reader, chunks);

    for result in query.lines() {
        let line = result?;

        if intersects(header, &line, &region)? {
            println!("{line}");
        }
    }

    Ok(())
}
