//! Prints an alignment file in the SAM format.
//!
//! Reference sequences in the FASTA format are only required for CRAM inputs that require them.
//!
//! The result matches the output of `samtools view --no-PG --with-header [--reference <fasta-src>]
//! <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufWriter},
};

use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_sam::{self as sam, AlignmentWriter};
use noodles_util::alignment;

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let fasta_src = args.next();

    let mut builder = File::open(src).map(alignment::Reader::builder)?;

    if let Some(fasta_src) = fasta_src {
        let repository = IndexedReader::builder()
            .open(fasta_src)
            .map(fasta::Repository::new)?;

        builder = builder.set_reference_sequence_repository(repository);
    }

    let mut reader = builder.build()?;
    let header = reader.read_header()?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = sam::Writer::new(BufWriter::new(handle));

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        writer.write_alignment_record(&header, record.as_ref())?;
    }

    Ok(())
}
