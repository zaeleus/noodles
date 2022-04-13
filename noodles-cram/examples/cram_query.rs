//! Queries a CRAM file with a given region.
//!
//! The input CRAM must have an associated CRAI in the same directory.
//!
//! The result matches the output of `samtools view [--reference <fasta-src>] <src> <region`.

use std::{
    env,
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use noodles_cram::{self as cram, crai};
use noodles_fasta::{self as fasta, fai, repository::adapters::IndexedReader};
use noodles_sam as sam;
use sam::AlignmentWriter;

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

    let src = args.next().map(PathBuf::from).expect("missing src");
    let region = args.next().expect("missing region").parse()?;
    let fasta_src = args.next();

    let repository = fasta_src
        .map(create_reference_sequence_repository)
        .transpose()?
        .unwrap_or_default();

    let mut reader = File::open(&src).map(cram::Reader::new)?;
    reader.read_file_definition()?;
    let header = reader.read_file_header()?.parse()?;

    let index = crai::read(src.with_extension("cram.crai"))?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = sam::Writer::new(BufWriter::new(handle));

    let query = reader.query(&repository, &header, &index, &region)?;

    for result in query {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    writer.into_inner().flush()?;

    Ok(())
}
