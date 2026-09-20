//! Queries an alignment file with the given region.
//!
//! The input must have an associated index in the same directory.
//!
//! The result matches the output of `samtools view [--reference <fasta-src>] <src> <region>`.

use std::{env, pin::Pin};

use futures::{Stream, TryStreamExt};
use noodles_fasta as fasta;
use noodles_sam as sam;
use noodles_util::alignment;
use tokio::io::{self, AsyncWriteExt};

const UNMAPPED: &str = "*";

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let raw_region = args.next().expect("missing region");
    let fasta_src = args.next();

    let mut builder = alignment::r#async::io::reader::Builder::default();

    if let Some(fasta_src) = fasta_src {
        let repository = fasta::io::indexed_reader::Builder::default()
            .build_from_path(fasta_src)
            .map(fasta::repository::adapters::IndexedReader::new)
            .map(fasta::Repository::new)?;

        builder = builder.set_reference_sequence_repository(repository);
    }

    let mut reader = builder.build_from_path(&src).await?;
    let header = reader.read_header().await?;

    let index = alignment::r#async::fs::read_associated_index(&src).await?;

    let mut query: Pin<Box<dyn Stream<Item = io::Result<_>>>> = if raw_region == UNMAPPED {
        reader.query_unmapped(&header, &index).await.map(Box::pin)?
    } else {
        let region = raw_region.parse()?;
        reader.query(&header, &index, &region).map(Box::pin)?
    };

    let mut writer = sam::r#async::io::Writer::new(io::stdout());

    while let Some(record) = query.try_next().await? {
        writer.write_alignment_record(&header, &record).await?;
    }

    writer.get_mut().shutdown().await?;

    Ok(())
}
