//! Builds a FASTA index.
//!
//! This writes the output to stdout rather than `<src>.fai`.
//!
//! The result matches the output of `samtools faidx <src>`.

use std::{env, io};

use noodles_fasta::{self as fasta, fai};

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let index = fasta::index(src)?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = fai::Writer::new(handle);

    writer.write_index(&index)?;

    Ok(())
}
