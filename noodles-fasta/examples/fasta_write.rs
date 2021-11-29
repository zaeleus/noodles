//! Creates a new FASTA file.
//!
//! This writes a single FASTA record to stdout.

use std::io;

use noodles_fasta as fasta;

fn main() -> io::Result<()> {
    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = fasta::Writer::new(handle);

    let definition = fasta::record::Definition::new("sq0", None);
    let sequence = b"ACGT".iter().cycle().take(256).copied().collect();
    let record = fasta::Record::new(definition, sequence);

    writer.write_record(&record)?;

    Ok(())
}
