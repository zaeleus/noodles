//! Queries a FASTA with a given reference sequence name.
//!
//! The input FASTA must have both gzi and fai in the same directory.
//!
//! The result is similar to the output of `samtools faidx --length 80 <src>
//! <reference-sequence-name>`.

use std::{
    env,
    fs::File,
    io,
    path::PathBuf,
};

use noodles_bgzf as bgzf;
use noodles_fasta::{self as fasta, fai};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut gzi_reader = File::open(src.with_extension("gz.gzi"))
        .map(bgzf::gzi::Reader::new)?;
    let gzi = gzi_reader.read_index()?;

    let mut reader = File::open(&src)
        .map(|f| bgzf::Reader::with_gzi(f, gzi))
        .map(fasta::Reader::new)?;

    let index = fai::read(src.with_extension("gz.fai"))?;
    let region = raw_region.parse()?;

    let record = reader.query(&index, &region)?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = fasta::Writer::new(handle);

    writer.write_record(&record)?;

    Ok(())
}
