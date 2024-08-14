//! Prints an alignment file in the SAM format.
//!
//! Reference sequences in the FASTA format are only required for CRAM inputs that require them.
//!
//! The result matches the output of `samtools view --no-PG --with-header [--reference <fasta-src>]
//! <src>`.

use futures::TryStreamExt;
use std::env;
use tokio::io::{self, BufWriter};

use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_sam::r#async::io::Writer;
use noodles_util::alignment::r#async::io::reader::Builder;

#[tokio::main]
async fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let fasta_src = args.next();

    let mut builder = Builder::default();

    if let Some(fasta_src) = fasta_src {
        let repository = fasta::io::indexed_reader::Builder::default()
            .build_from_path(fasta_src)
            .map(IndexedReader::new)
            .map(fasta::Repository::new)?;

        builder = builder.set_reference_sequence_repository(repository);
    }

    let mut reader = if src == "-" {
        let stdin = io::stdin();
        builder.build_from_reader(stdin).await?
    } else {
        builder.build_from_path(src).await?
    };

    let header = reader.read_header().await?;

    let stdout = io::stdout();
    let mut writer = Writer::new(BufWriter::new(stdout));

    writer.write_header(&header).await?;

    while let Some(record) = reader.records(&header).try_next().await? {
        writer.write_alignment_record(&header, &record).await?;
    }

    Ok(())
}
