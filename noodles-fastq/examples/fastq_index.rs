//! Builds and writes a FASTQ index.
//!
//! This writes the output to stdout rather than `<src>.fai`.
//!
//! The result matches the output of `samtools fqidx <src>`.

use std::{env, io};

use noodles_fastq::{self as fastq, fai};

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let index = fastq::io::index(src)?;

    let stdout = io::stdout().lock();
    let mut writer = fai::io::Writer::new(stdout);

    for record in &index {
        writer.write_record(record)?;
    }

    Ok(())
}
