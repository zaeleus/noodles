//! Replaces the SAM header of a BAM file.
//!
//! This is similar to the functionality of `samtools reheader`.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::{env, fs::File, io};

use noodles_bam as bam;
use noodles_sam::header::record::value::{
    Map,
    map::{Program, builder::BuildError, program::tag},
};

fn build_self_program() -> Result<Map<Program>, BuildError> {
    let args = env::args().collect::<Vec<_>>().join(" ");

    Map::builder()
        .insert(tag::NAME, env!("CARGO_BIN_NAME"))
        .insert(tag::VERSION, env!("CARGO_PKG_VERSION"))
        .insert(tag::COMMAND_LINE, args)
        .build()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let mut header = reader.read_header()?;

    let pg = build_self_program()?;
    header.programs_mut().add("noodles", pg)?;

    header.add_comment("a comment added by noodles-bam");

    let stdout = io::stdout().lock();
    let mut writer = bam::io::Writer::new(stdout);

    writer.write_header(&header)?;

    io::copy(reader.get_mut(), writer.get_mut())?;

    Ok(())
}
