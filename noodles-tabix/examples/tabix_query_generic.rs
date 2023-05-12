//! Queries an indexed TSV with the given region.
//!
//! The input must have an associated tabix index in the same directory.
//!
//! The result matches the output of `tabix <src> <region>`.

use std::{env, fs::File, io};

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi as csi;
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
    let records = csi::io::FilteredLines::new(query, header, region);

    for result in records {
        let line = result?;
        println!("{line}");
    }

    Ok(())
}
