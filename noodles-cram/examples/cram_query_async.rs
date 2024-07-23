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
use tokio::{fs::File, io};

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

    let mut reader = File::open(&src).await.map(cram::AsyncReader::new)?;
    let header = reader.read_header().await?;

    let index = crai::r#async::read(src.with_extension("cram.crai")).await?;
    let mut query = reader.query(&repository, &header, &index, &region)?;

    let mut writer = sam::AsyncWriter::new(io::stdout());

    while let Some(record) = query.try_next().await? {
        let r = record.try_into_alignment_record(&header)?;
        writer.write_alignment_record(&header, &r).await?;
    }

    Ok(())
}
