//! Rewrites an alignment format to another alignment format asynchronously.
//!
//! The output format is determined from the extension of the destination.

use std::env;

use futures::TryStreamExt;
use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_util::alignment;
use tokio::io;

#[tokio::main]
async fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let dst = args.next().expect("missing dst");
    let fasta_src = args.next();

    let repository = fasta_src
        .map(|src| fasta::io::indexed_reader::Builder::default().build_from_path(src))
        .transpose()?
        .map(IndexedReader::new)
        .map(fasta::Repository::new)
        .unwrap_or_default();

    let mut reader = alignment::r#async::io::reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .build_from_path(src)
        .await?;

    let header = reader.read_header().await?;

    let mut writer = alignment::r#async::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(dst)
        .await?;

    writer.write_header(&header).await?;

    while let Some(record) = reader.records(&header).try_next().await? {
        writer.write_record(&header, &record).await?;
    }

    writer.finish(&header).await?;

    Ok(())
}
