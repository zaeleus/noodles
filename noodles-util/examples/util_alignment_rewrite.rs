//! Rewrites an alignment format to another alignment format.
//!
//! The output format is determined from the extension of the destination.

use std::{env, io, path::Path};

use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_util::alignment::{self, Format};

fn detect_format_from_extension<P>(path: P) -> Option<Format>
where
    P: AsRef<Path>,
{
    match path.as_ref().extension().and_then(|ext| ext.to_str()) {
        Some("sam") => Some(Format::Sam),
        Some("bam") => Some(Format::Bam),
        Some("cram") => Some(Format::Cram),
        _ => None,
    }
}

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

    let mut reader = alignment::reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .build_from_path(src)?;

    let header = reader.read_header()?;

    let output_format = detect_format_from_extension(&dst).expect("invalid dst extension");

    let mut writer = alignment::writer::Builder::default()
        .set_format(output_format)
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
