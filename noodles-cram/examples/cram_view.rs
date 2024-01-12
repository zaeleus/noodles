//! Prints a CRAM file in the SAM format.
//!
//! Reference sequences in the FASTA format are only required for inputs that require them.
//!
//! The result matches the output of `samtools view [--reference <fasta-src>] <src>`.

use std::{
    env,
    io::{self, BufWriter},
};

use noodles_cram as cram;
use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_sam::{self as sam, alignment::io::Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let fasta_src = args.next();

    let reference_sequence_repository = fasta_src
        .map(|src| fasta::indexed_reader::Builder::default().build_from_path(src))
        .transpose()?
        .map(IndexedReader::new)
        .map(fasta::Repository::new)
        .unwrap_or_default();

    let mut reader = cram::reader::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .build_from_path(src)?;

    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = sam::io::Writer::new(BufWriter::new(stdout));

    for result in reader.records(&header) {
        let record = result.and_then(|record| record.try_into_alignment_record(&header))?;
        writer.write_alignment_record(&header, &record)?;
    }

    Ok(())
}
