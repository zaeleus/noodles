//! Prints an alignment file in the SAM format.
//!
//! Reference sequences in the FASTA format are only required for CRAM inputs that require them.
//!
//! The result matches the output of `samtools view --no-PG --with-header [--reference <fasta-src>]
//! <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufWriter, Read},
};

use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_sam::{self as sam, alignment::io::Write};
use noodles_util::alignment;

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let fasta_src = args.next();

    let mut builder = alignment::io::reader::Builder::default();

    if let Some(fasta_src) = fasta_src {
        let repository = fasta::io::indexed_reader::Builder::default()
            .build_from_path(fasta_src)
            .map(IndexedReader::new)
            .map(fasta::Repository::new)?;

        builder = builder.set_reference_sequence_repository(repository);
    }

    let source: Box<dyn Read> = if src == "-" {
        Box::new(io::stdin().lock())
    } else {
        File::open(src).map(Box::new)?
    };

    let mut reader = builder.build_from_reader(source)?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = sam::io::Writer::new(BufWriter::new(stdout));

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    Ok(())
}
