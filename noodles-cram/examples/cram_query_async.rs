//! Queries a CRAM file with a given region.
//!
//! The input CRAM must have an associated CRAI in the same directory.
//!
//! The result matches the output of `samtools view [--reference <fasta-src>] <src> <region>`.

use std::{env, path::PathBuf, pin::Pin};

use futures::{Stream, TryStreamExt};
use noodles_cram::{self as cram, crai};
use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_sam as sam;
use tokio::io;

const UNMAPPED: &str = "*";

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().map(PathBuf::from).expect("missing src");
    let raw_region = args.next().expect("missing region");
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

    let mut records: Pin<Box<dyn Stream<Item = io::Result<sam::alignment::RecordBuf>>>> =
        if raw_region == UNMAPPED {
            reader.query_unmapped(&header, &index).await.map(Box::pin)?
        } else {
            let region = raw_region.parse()?;
            Box::pin(reader.query(&header, &index, &region)?)
        };

    let mut writer = sam::r#async::io::Writer::new(io::stdout());

    while let Some(record) = records.try_next().await? {
        writer.write_alignment_record(&header, &record).await?;
    }

    Ok(())
}
