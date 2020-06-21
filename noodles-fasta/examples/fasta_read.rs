//! Prints the reference sequence names and lengths of all the records in a FASTA file.

use std::{env, fs::File, io::BufReader};

use noodles_fasta as fasta;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(fasta::Reader::new)?;

    let mut definition_buf = String::new();
    let mut sequence_buf = Vec::new();

    loop {
        definition_buf.clear();

        match reader.read_definition(&mut definition_buf) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => return Err(e.into()),
        }

        let definition: fasta::record::Definition = definition_buf.parse()?;

        sequence_buf.clear();
        reader.read_sequence(&mut sequence_buf)?;

        println!(
            "{}\t{}",
            definition.reference_sequence_name(),
            sequence_buf.len()
        );
    }

    Ok(())
}
