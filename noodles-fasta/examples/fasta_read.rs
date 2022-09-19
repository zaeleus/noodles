//! Prints the reference sequence names and lengths of all the records in a FASTA file.

use std::{env, io};

use noodles_fasta as fasta;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = fasta::reader::Builder::default().build_from_path(src)?;

    for result in reader.records() {
        let record = result?;

        println!("{}\t{}", record.name(), record.sequence().len());
    }

    Ok(())
}
