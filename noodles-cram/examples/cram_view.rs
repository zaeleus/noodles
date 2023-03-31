//! Prints a CRAM file in the SAM format.
//!
//! Reference sequences in the FASTA format are only required for inputs that require them.
//!
//! The result matches the output of `samtools view [--reference <fasta-src>] <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufWriter},
};

use noodles_cram as cram;
use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_sam::{self as sam, AlignmentWriter};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let fasta_src = args.next();

    let repository = fasta_src
        .map(|src| fasta::indexed_reader::Builder::default().build_from_path(src))
        .transpose()?
        .map(IndexedReader::new)
        .map(fasta::Repository::new)
        .unwrap_or_default();

    let mut reader = File::open(src).map(cram::Reader::new)?;
    reader.read_file_definition()?;

    let header = reader.read_file_header()?;

    let stdout = io::stdout().lock();
    let mut writer = sam::Writer::new(BufWriter::new(stdout));

    for result in reader.records(&repository, &header) {
        let record = result.and_then(|record| record.try_into_alignment_record(&header))?;
        writer.write_alignment_record(&header, &record)?;
    }

    Ok(())
}
