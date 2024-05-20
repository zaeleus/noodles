use std::env;

use futures::TryStreamExt;
use noodles_cram as cram;
use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_sam as sam;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let fasta_src = args.next();

    let repository = fasta_src
        .map(|src| fasta::io::indexed_reader::Builder::default().build_from_path(src))
        .transpose()?
        .map(IndexedReader::new)
        .map(fasta::Repository::new)
        .unwrap_or_default();

    let mut reader = File::open(src).await.map(cram::AsyncReader::new)?;
    reader.read_file_definition().await?;

    let header = reader.read_file_header().await?.parse()?;

    let mut records = reader.records(&repository, &header);

    let mut writer = sam::AsyncWriter::new(io::stdout());

    while let Some(record) = records.try_next().await? {
        let record = record.try_into_alignment_record(&header)?;
        writer.write_alignment_record(&header, &record).await?;
    }

    Ok(())
}
