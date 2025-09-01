//! Queries a bgzipped BED file with a given region.
//!
//! The input bgzipped BED file must have an associated tabix index (TBI) in the same directory.

use std::{env, fs::File, io};

use noodles_bed as bed;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::{self as csi, BinningIndex};
use noodles_tabix as tabix;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let region: Region = args.next().expect("missing region").parse()?;

    let index_src = format!("{src}.tbi");
    let index = tabix::fs::read(index_src)?;

    let header = index.header().expect("missing tabix header");
    let reference_sequence_id = header
        .reference_sequence_names()
        .get_index_of(region.name())
        .expect("invalid reference sequence name");

    let mut decoder = File::open(src).map(bgzf::io::Reader::new)?;
    let chunks = index.query(reference_sequence_id, region.interval())?;
    let query = csi::io::Query::new(&mut decoder, chunks);

    let mut reader = bed::io::Reader::<3, _>::new(query);
    let mut record = bed::Record::default();

    let stdout = io::stdout().lock();
    let mut writer = bed::io::Writer::<3, _>::new(stdout);

    while reader.read_record(&mut record)? != 0 {
        let start = record.feature_start()?;
        let end = record.feature_end().expect("missing feature end")?;
        let interval = (start..=end).into();

        if !region.interval().intersects(interval) {
            continue;
        }

        writer.write_record(&record)?;
    }

    Ok(())
}
