//! Add a PG tag to the SAM header of a BAM file.
//!
//! This is a demonstration of how to construct a full PG entry.
//!
//! Verify the output by piping to `samtools view --with-header`.

use std::{env, io};

use noodles_bam as bam;
use noodles_sam::header::record::value::{
    map::{Program, Tag},
    Map,
};
use noodles_sam::Header;

fn add_pg(mut header: Header) -> Header {
    const NAME: &str = "bam_add_tag_to_header";
    const VERSION: &str = env!("CARGO_PKG_VERSION");

    // standard PG tags PN, VN, PP, CL
    let pn = match Tag::try_from([b'P', b'N']) {
        Ok(Tag::Other(tag)) => tag,
        _ => unreachable!(),
    };
    let vn = match Tag::try_from([b'V', b'N']) {
        Ok(Tag::Other(tag)) => tag,
        _ => unreachable!(),
    };
    let pp = match Tag::try_from([b'P', b'P']) {
        Ok(Tag::Other(tag)) => tag,
        _ => unreachable!(),
    };
    let cl = match Tag::try_from([b'C', b'L']) {
        Ok(Tag::Other(tag)) => tag,
        _ => unreachable!(),
    };
    // the command-line to insert into the CL tag
    let cmd_str = env::args().collect::<Vec<_>>().join(" ");

    let program = Map::<Program>::builder().insert(pn, NAME);

    let program = if let Some(last_pg) = header.programs().keys().last() {
        program.insert(pp, last_pg.clone())
    } else {
        program
    };

    let program = program
        .insert(vn, VERSION)
        .insert(cl, cmd_str)
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
