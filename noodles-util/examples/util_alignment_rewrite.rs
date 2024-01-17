//! Rewrites an alignment format to another alignment format.
//!
//! The output format is determined from the extension of the destination.

use std::{env, io};

use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_util::alignment;

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let dst = args.next().expect("missing dst");
    let fasta_src = args.next();

    let repository = fasta_src
        .map(|src| fasta::indexed_reader::Builder::default().build_from_path(src))
        .transpose()?
        .map(IndexedReader::new)
        .map(fasta::Repository::new)
        .unwrap_or_default();

    let mut reader = alignment::io::reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .build_from_path(src)?;

    let header = reader.read_header()?;

    let mut writer = alignment::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(dst)?;

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    writer.finish(&header)?;

    Ok(())
}
