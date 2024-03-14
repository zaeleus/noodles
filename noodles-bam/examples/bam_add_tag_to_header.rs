//! Add a PG tag to the SAM header of a BAM file.
//!
//! This is a demonstration of how to construct a full PG entry.
//!
//! Verify the output by piping to `samtools view --with-header`.

use std::{env, io};

use noodles_bam as bam;
use noodles_sam::header::record::value::{
    map::{program::tag, Program},
    Map,
};
use noodles_sam::Header;

fn add_pg(mut header: Header) -> Header {
    const NAME: &str = "bam_add_tag_to_header";
    const VERSION: &str = env!("CARGO_PKG_VERSION");

    // the command-line to insert into the CL tag
    let cmd_str = env::args().collect::<Vec<_>>().join(" ");

    let program = Map::<Program>::builder().insert(tag::NAME, NAME);

    let program = if let Some(last_pg) = header.programs().keys().last() {
        program.insert(tag::PREVIOUS_PROGRAM_ID, last_pg.clone())
    } else {
        program
    };

    let program = program
        .insert(tag::VERSION, VERSION)
        .insert(tag::COMMAND_LINE, cmd_str)
        .build()
        .unwrap();
    header.programs_mut().insert(NAME.into(), program);

    header
}

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bam::io::reader::Builder.build_from_path(src)?;

    let header = add_pg(reader.read_header()?);

    let stdout = io::stdout().lock();
    let mut writer = bam::io::Writer::new(stdout);

    writer.write_header(&header)?;

    io::copy(reader.get_mut(), writer.get_mut())?;

    Ok(())
}
