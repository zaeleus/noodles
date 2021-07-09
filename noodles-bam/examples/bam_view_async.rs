//! Prints a BAM file in the SAM format.
//!
//! The result matches the output of `samtools view <src>`.

use std::env;

use noodles_bam as bam;
use noodles_sam as sam;
use tokio::fs::File;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(bam::AsyncReader::new)?;
    let header: sam::Header = reader.read_header().await?.parse()?;
    reader.read_reference_sequences().await?;

    let mut record = bam::Record::default();

    loop {
        match reader.read_record(&mut record).await {
            Ok(0) => break,
            Ok(_) => {
                let sam_record = record.try_into_sam_record(header.reference_sequences())?;
                println!("{}", sam_record);
            }
            Err(e) => return Err(e.into()),
        }
    }

    Ok(())
}
