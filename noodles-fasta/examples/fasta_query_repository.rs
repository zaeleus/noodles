//! Queries a FASTA with a given region.
//!
//! The input FASTA must have an index in the same directory. It uses a repository cache that
//! allows for more efficient (repeated) lookups of a sequence.
//!
//! The result is similar to the output of `samtools faidx --length 80 <src> <region>`.

use std::{env, io};

use noodles_fasta::{self as fasta, record::Definition, repository::adapters::IndexedReader};

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let name = args.next().expect("missing name");

    let reader = fasta::io::indexed_reader::Builder::default().build_from_path(src)?;
    let adapter = IndexedReader::new(reader);
    let repository = fasta::Repository::new(adapter);

    let sequence = repository
        .get(name.as_bytes())
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid name"))?;

    let stdout = io::stdout().lock();
    let mut writer = fasta::io::Writer::new(stdout);

    let record = fasta::Record::new(Definition::new(name, None), sequence);
    writer.write_record(&record)?;

    Ok(())
}
