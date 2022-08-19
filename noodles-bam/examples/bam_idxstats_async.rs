//! Prints BAM index statistics.
//!
//! The data is read from both the SAM header reference sequences and BAM index. It is printed as a
//! tab-delimited record with the following columns for a region: reference sequence name,
//! reference sequence length, number of mapped records, and number of unmapped records.
//!
//! The result matches the output of `samtools idxstats <src>`.

use std::{env, path::PathBuf};

use noodles_bam::{self as bam, bai};
use noodles_csi::{binning_index::ReferenceSequenceExt, BinningIndex};
use noodles_sam as sam;
use tokio::fs::File;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).map(PathBuf::from).expect("missing src");

    let mut reader = File::open(&src).await.map(bam::AsyncReader::new)?;
    let header: sam::Header = reader.read_header().await?.parse()?;

    let index = bai::r#async::read(src.with_extension("bam.bai")).await?;

    for (reference_sequence, index_reference_sequence) in header
        .reference_sequences()
        .values()
        .zip(index.reference_sequences())
    {
        let (mapped_record_count, unmapped_record_count) = index_reference_sequence
            .metadata()
            .map(|m| (m.mapped_record_count(), m.unmapped_record_count()))
            .unwrap_or_default();

        println!(
            "{}\t{}\t{}\t{}",
            reference_sequence.name(),
            reference_sequence.length(),
            mapped_record_count,
            unmapped_record_count
        );
    }

    let unmapped_record_count = index.unplaced_unmapped_record_count().unwrap_or_default();
    println!("*\t0\t0\t{}", unmapped_record_count);

    Ok(())
}
