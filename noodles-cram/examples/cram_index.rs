//! Builds and writes a CRAM index from a CRAM file.
//!
//! This writes the output to stdout rather than `<src>.crai`.
//!
//! The output is similar to the output of `samtools index <src>`.

use std::{env, io};

use noodles_cram::{self as cram, crai};

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let index = cram::fs::index(src)?;

    let stdout = io::stdout().lock();
    let mut writer = crai::io::Writer::new(stdout);

    writer.write_index(&index)?;

    Ok(())
}
