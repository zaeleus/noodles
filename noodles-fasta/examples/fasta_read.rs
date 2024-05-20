//! Prints the reference sequence names and lengths of all the records in a FASTA file.

use std::{env, io, str};

use noodles_fasta as fasta;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = fasta::io::reader::Builder.build_from_path(src)?;

    for result in reader.records() {
        let record = result?;

        let name = str::from_utf8(record.name())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let length = record.sequence().len();

        println!("{name}\t{length}");
    }

    Ok(())
}
