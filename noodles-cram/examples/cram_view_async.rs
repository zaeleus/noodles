use std::env;

use futures::TryStreamExt;
use noodles_cram as cram;
use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use tokio::fs::File;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let fasta_src = args.next();

    let repository = fasta_src
        .map(|src| IndexedReader::builder().open(src))
        .transpose()?
        .map(fasta::Repository::new)
        .unwrap_or_default();

    let mut reader = File::open(src).await.map(cram::AsyncReader::new)?;
    reader.read_file_definition().await?;

    let header = reader.read_file_header().await?.parse()?;

    let mut records = reader.records(&repository, &header);

    while let Some(record) = records.try_next().await? {
        let sam_record = record.try_into_sam_record(&header)?;
        println!("{}", sam_record);
    }

    Ok(())
}
