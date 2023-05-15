//! Creates a new FASTQ file.
//!
//! This writes a single FASTQ record to stdout.

use noodles_fastq::{self as fastq, record::Definition};
use tokio::io;

#[tokio::main]
async fn main() -> io::Result<()> {
    let mut writer = fastq::AsyncWriter::new(io::stdout());

    let record = fastq::Record::new(Definition::new("r0", ""), "ACGT", "NDLS");
    writer.write_record(&record).await?;

    Ok(())
}
