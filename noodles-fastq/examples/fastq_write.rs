//! Creates a new FASTQ file.
//!
//! This writes a single FASTQ record to stdout.

use std::io;

use noodles_fastq as fastq;

fn main() -> io::Result<()> {
    let stdout = io::stdout().lock();
    let mut writer = fastq::Writer::new(stdout);

    let record = fastq::Record::new("r0", "ACGT", "NDLS");
    writer.write_record(&record)?;

    Ok(())
}
