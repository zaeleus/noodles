//! Queries a CRAM file with a given region.
//!
//! The input CRAM must have an associated CRAI in the same directory.
//!
//! The result matches the output of `samtools view [--reference <fasta-src>] <src> <region`.

use std::{
    env,
    fs::File,
    io::{self, BufWriter, Write},
    path::PathBuf,
};

use noodles_cram::{self as cram, crai};
use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_sam as sam;
use sam::AlignmentWriter;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().map(PathBuf::from).expect("missing src");
    let region = args.next().expect("missing region").parse()?;
    let fasta_src = args.next();

    let repository = fasta_src
        .map(|src| fasta::indexed_reader::Builder::default().build_from_path(src))
        .transpose()?
        .map(IndexedReader::new)
        .map(fasta::Repository::new)
        .unwrap_or_default();

    let mut reader = File::open(&src).map(cram::Reader::new)?;
    reader.read_file_definition()?;
    let header = reader.read_file_header()?.parse()?;

    let index = crai::read(src.with_extension("cram.crai"))?;

    let stdout = io::stdout().lock();
    let mut writer = sam::Writer::new(BufWriter::new(stdout));

    let query = reader.query(&repository, &header, &index, &region)?;

    for result in query {
        let record = result.and_then(|record| record.try_into_alignment_record(&header))?;
        writer.write_alignment_record(&header, &record)?;
    }

    writer.into_inner().flush()?;

    Ok(())
}
