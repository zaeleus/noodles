//! Queries a FASTA with a given region.
//!
//! The input FASTA must have an index in the same directory.
//!
//! The result is similar to the output of `samtools faidx --length 80 <src> <region>`.

use std::{env, io};

use noodles_fasta as fasta;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = fasta::indexed_reader::Builder::default().build_from_path(src)?;

    let region = raw_region.parse()?;
    let record = reader.query(&region)?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = fasta::Writer::new(handle);

    writer.write_record(&record)?;

    Ok(())
}
