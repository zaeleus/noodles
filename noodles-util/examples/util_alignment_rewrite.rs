//! Rewrites an alignment format to another alignment format.
//!
//! The output format is determined from the extension of the destination.

use std::{
    env,
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, BufReader, BufWriter},
    path::{Path, PathBuf},
};

use noodles_fasta::{self as fasta, fai, repository::adapters::IndexedReader};
use noodles_util::alignment::{self, Format};

fn detect_format_from_extension<P>(path: P) -> Option<Format>
where
    P: AsRef<Path>,
{
    match path.as_ref().extension().and_then(|ext| ext.to_str()) {
        Some("sam") => Some(Format::Sam),
        Some("bam") => Some(Format::Bam),
        Some("cram") => Some(Format::Cram),
        _ => None,
    }
}

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
    let dst = args.next().expect("missing dst");
    let fasta_src = args.next();

    let repository = fasta_src
        .map(create_reference_sequence_repository)
        .transpose()?
        .unwrap_or_default();

    let mut reader = File::open(src).and_then(|f| {
        alignment::Reader::builder(f)
            .set_reference_sequence_repository(repository.clone())
            .build()
    })?;

    let header = reader.read_header()?;

    let output_format = detect_format_from_extension(&dst).expect("invalid dst extension");

    let mut writer = File::create(dst).map(|f| {
        alignment::Writer::builder(BufWriter::new(f))
            .set_format(output_format)
            .set_reference_sequence_repository(repository)
            .build()
    })?;

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    writer.finish(&header)?;

    Ok(())
}
