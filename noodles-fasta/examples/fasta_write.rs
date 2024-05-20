//! Creates a new FASTA file.
//!
//! This writes a single FASTA record to stdout.

use std::io;

use noodles_fasta as fasta;

fn main() -> io::Result<()> {
    let stdout = io::stdout().lock();
    let mut writer = fasta::io::Writer::new(stdout);

    let definition = fasta::record::Definition::new("sq0", None);
    let sequence = fasta::record::Sequence::from(b"ACGT".repeat(64));
    let record = fasta::Record::new(definition, sequence);

    writer.write_record(&record)?;

    Ok(())
}
