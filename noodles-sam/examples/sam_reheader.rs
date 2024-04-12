//! Replaces the header of a SAM file.
//!
//! This is similar to the functionality of `samtools reheader`.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::{env, io};

use noodles_sam::{
    self as sam,
    header::record::value::{
        map::{builder::BuildError, program::tag, Program},
        Map,
    },
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

    let mut reader = sam::io::reader::Builder::default().build_from_path(src)?;
    let mut header = reader.read_header()?;

    let pg = build_self_program()?;
    header.programs_mut().add("noodles", pg)?;

    header.add_comment("a comment added by noodles-sam");

    let stdout = io::stdout().lock();
    let mut writer = sam::io::Writer::new(stdout);

    writer.write_header(&header)?;

    io::copy(reader.get_mut(), writer.get_mut())?;

    Ok(())
}
