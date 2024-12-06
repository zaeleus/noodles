//! Calculates the read depth of each position in a region.
//!
//! The results match the output of `samtools depth -r <region> <src>`.

use std::{
    env,
    io::{self, BufWriter, Write},
};

use noodles_core::Region;
use noodles_fasta as fasta;
use noodles_util::alignment::{self, iter::Depth};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let region: Region = args.next().expect("missing region").parse()?;
    let fasta_src = args.next();

    let mut builder = alignment::io::indexed_reader::Builder::default();

    if let Some(fasta_src) = fasta_src {
        let repository = fasta::io::indexed_reader::Builder::default()
            .build_from_path(fasta_src)
            .map(fasta::repository::adapters::IndexedReader::new)
            .map(fasta::Repository::new)?;

        builder = builder.set_reference_sequence_repository(repository);
    }

    let mut reader = builder.build_from_path(src)?;
    let header = reader.read_header()?;

    let query = reader.query(&header, &region)?;
    let pileup = Depth::new(&header, query);

    let stdout = io::stdout().lock();
    let mut writer = BufWriter::new(stdout);

    let reference_sequence_name = region.name();

    for result in pileup {
        let (position, depth) = result?;

        if !region.interval().contains(position) {
            continue;
        }

        writeln!(writer, "{reference_sequence_name}\t{position}\t{depth}")?;
    }

    Ok(())
}
