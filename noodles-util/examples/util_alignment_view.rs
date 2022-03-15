//! Prints an alignment file in the SAM format.
//!
//! Reference sequences in the FASTA format are only required for CRAM inputs that require them.
//!
//! The result matches the output of `samtools view --no-PG --with-header [--reference <fasta-src>]
//! <src>`.

use std::{
    env,
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, BufReader},
    path::{Path, PathBuf},
};

use fasta::repository::adapters::IndexedReader;
use noodles_fasta::{self as fasta, fai};
use noodles_sam::{self as sam, AlignmentWriter};
use noodles_util::alignment;

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

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let fasta_src = args.next();

    let mut builder = File::open(src).map(alignment::Reader::builder)?;

    if let Some(fasta_src) = fasta_src {
        let repository = create_reference_sequence_repository(fasta_src)?;
        builder = builder.set_reference_sequence_repository(repository);
    }

    let mut reader = builder.build()?;
    let header = reader.read_header()?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = sam::Writer::new(handle);

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    Ok(())
}
