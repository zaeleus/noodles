//! Queries a FASTA with a given reference sequence name.
//!
//! The input FASTA must have an index in the same directory.
//!
//! The result is similar to the output of `samtools faidx --length 80 <src>
//! <reference-sequence-name>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
    path::PathBuf,
};

use noodles_fasta::{self as fasta, fai};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = File::open(&src)
        .map(BufReader::new)
        .map(fasta::Reader::new)?;

    let index = fai::read(src.with_extension("fa.fai"))?;
    let region = raw_region.parse()?;

    let record = reader.query(&index, &region)?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = fasta::Writer::new(handle);

    writer.write_record(&record)?;

    Ok(())
}
