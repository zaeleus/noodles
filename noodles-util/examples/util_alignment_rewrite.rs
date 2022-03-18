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

    let mut reader = File::open(src).and_then(|f| alignment::Reader::builder(f).build())?;
    let header = reader.read_header()?;

    let mut builder = File::create(&dst)
        .map(BufWriter::new)
        .map(alignment::Writer::builder)?;

    let output_format = detect_format_from_extension(dst).expect("invalid dst extension");
    builder = builder.set_format(output_format);

    if let Some(fasta_src) = fasta_src {
        let repository = create_reference_sequence_repository(fasta_src)?;
        builder = builder.set_reference_sequence_repository(repository);
    }

    let mut writer = builder.build();

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
