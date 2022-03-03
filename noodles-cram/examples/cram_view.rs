//! Prints a CRAM file in the SAM format.
//!
//! The result matches the output of `samtools view <src>`.

use std::{
    env,
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, BufReader},
    path::{Path, PathBuf},
};

use noodles_cram as cram;
use noodles_fasta::{self as fasta, fai, repository::adapters::IndexedReader};
use noodles_sam as sam;

fn push_ext<S>(path: PathBuf, ext: S) -> PathBuf
where
    S: AsRef<OsStr>,
{
    let mut s = OsString::from(path);
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

fn create_reference_sequence_repository<P>(src: P) -> io::Result<fasta::Repository>
where
    P: AsRef<Path>,
{
    let src = src.as_ref();

    let reader = File::open(src)
        .map(BufReader::new)
        .map(fasta::Reader::new)?;

    let fai_src = push_ext(src.to_path_buf(), "fai");
    let index = fai::read(fai_src)?;

    let adapter = IndexedReader::new(reader, index);

    Ok(fasta::Repository::new(adapter))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let fasta_src = args.next().expect("missing fasta_src");
    let src = args.next().expect("missing src");

    let reference_sequence_repository = create_reference_sequence_repository(fasta_src)?;

    let mut reader = File::open(src).map(cram::Reader::new)?;
    reader.read_file_definition()?;

    let header: sam::Header = reader.read_file_header()?.parse()?;

    for result in reader.records(&reference_sequence_repository, &header) {
        let record = result?;
        let sam_record = record.try_into_sam_record(header.reference_sequences())?;
        println!("{}", sam_record);
    }

    Ok(())
}
