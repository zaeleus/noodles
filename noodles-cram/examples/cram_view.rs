//! Prints a CRAM file in the SAM format.
//!
//! Reference sequences in the FASTA format are only required for inputs that require them.
//!
//! The result matches the output of `samtools view [--reference <fasta-src>] <src>`.

use std::{
    env,
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, BufReader, BufWriter},
    path::{Path, PathBuf},
};

use noodles_cram as cram;
use noodles_fasta::{self as fasta, fai, repository::adapters::IndexedReader};
use noodles_sam::{self as sam, AlignmentWriter};

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

    let src = args.next().expect("missing src");
    let fasta_src = args.next();

    let repository = if let Some(src) = fasta_src {
        create_reference_sequence_repository(src)?
    } else {
        fasta::Repository::default()
    };

    let mut reader = File::open(src).map(cram::Reader::new)?;
    reader.read_file_definition()?;

    let header: sam::Header = reader.read_file_header()?.parse()?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = sam::Writer::new(BufWriter::new(handle));

    for result in reader.records(&repository, &header) {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    Ok(())
}
