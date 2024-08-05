//! Creates a new BED3 file.
//!
//! This writes a single BED3 record to stdout.

use std::io;

use noodles_bed::{self as bed, feature::RecordBuf};
use noodles_core::Position;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let stdout = io::stdout().lock();
    let mut writer = bed::io::Writer::<3, _>::new(stdout);

    let record = RecordBuf::<3>::builder()
        .set_reference_sequence_name("sq0")
        .set_feature_start(Position::try_from(8)?)
        .set_feature_end(Position::try_from(13)?)
        .build();

    writer.write_feature_record(&record)?;

    Ok(())
}
