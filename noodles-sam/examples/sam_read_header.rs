//! Prints the header of a SAM file.
//!
//! The result is similar to the output of `samtools head <src>`.

use std::{env, io};

use noodles_sam as sam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = sam::io::reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = sam::io::Writer::new(stdout);
    writer.write_header(&header)?;

    Ok(())
}
