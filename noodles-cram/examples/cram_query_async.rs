//! Queries a CRAM file with a given region.
//!
//! The input CRAM must have an associated CRAI in the same directory.
//!
//! The result matches the output of `samtools view [--reference <fasta-src>] <src> <region>`.

use std::{env, path::PathBuf};

use futures::TryStreamExt;
use noodles_cram::{self as cram, crai};
use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_sam as sam;
use tokio::io;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().map(PathBuf::from).expect("missing src");
    let region = args.next().expect("missing region").parse()?;
    let fasta_src = args.next();

    let repository = fasta_src
        .map(|src| fasta::io::indexed_reader::Builder::default().build_from_path(src))
        .transpose()?
        .map(IndexedReader::new)
        .map(fasta::Repository::new)
        .unwrap_or_default();

    let mut reader = cram::r#async::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(&src)
        .await?;

    let header = reader.read_header().await?;

    let index = crai::r#async::read(src.with_extension("cram.crai")).await?;
    let mut query = reader.query(&header, &index, &region)?;

    let mut writer = sam::r#async::io::Writer::new(io::stdout());

    while let Some(record) = query.try_next().await? {
        writer.write_alignment_record(&header, &record).await?;
    }

    Ok(())
}
