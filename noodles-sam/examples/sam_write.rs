//! Creates a new SAM file.
//!
//! This writes a SAM header and one unmapped record to stdout.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::io;

use noodles_sam::{
    self as sam,
    alignment::{RecordBuf, io::Write},
    header::record::value::{Map, map::Program},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let stdout = io::stdout().lock();
    let mut writer = sam::io::Writer::new(stdout);

    let header = sam::Header::builder()
        .set_header(Default::default())
        .add_program("noodles-sam", Map::<Program>::default())
        .add_comment("an example SAM written by noodles-sam")
        .build();

    writer.write_header(&header)?;

    let record = RecordBuf::default();
    writer.write_alignment_record(&header, &record)?;

    Ok(())
}
