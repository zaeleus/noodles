//! Writes a BED-like file with a tabix index.
//!
//! This outputs to `out.tsv.gz` and `out.tsv.gz.tbi`. Use `tabix` or `tabix_query_generic` to test
//! querying the output.

use std::{fs::File, io::Write};

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_tabix as tabix;

const DST: &str = "out.tsv.gz";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let records: [(&str, Position, Position); 3] = [
        ("sq0", Position::try_from(8)?, Position::try_from(21)?),
        ("sq0", Position::try_from(13)?, Position::try_from(34)?),
        ("sq1", Position::try_from(5)?, Position::try_from(8)?),
    ];

    let mut writer = File::create(DST).map(bgzf::io::Writer::new)?;

    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::binning_index::index::header::Builder::bed().build());

    let mut start_position = writer.virtual_position();

    for (reference_sequence_name, start, end) in records {
        writeln!(writer, "{}\t{}\t{}", reference_sequence_name, start, end)?;

        let end_position = writer.virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        indexer.add_record(reference_sequence_name, start, end, chunk)?;

        start_position = end_position;
    }

    writer.finish()?;

    let index = indexer.build();

    let index_dst = format!("{DST}.tbi");
    let mut writer = File::create(index_dst).map(tabix::io::Writer::new)?;
    writer.write_index(&index)?;

    Ok(())
}
