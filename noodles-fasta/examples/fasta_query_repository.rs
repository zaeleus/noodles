use std::{
    env,
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, BufReader},
    path::PathBuf,
};

use noodles_fasta::{self as fasta, fai, record::Definition, repository::adapters::IndexedReader};

fn push_ext<S>(path: PathBuf, ext: S) -> PathBuf
where
    S: AsRef<OsStr>,
{
    let mut s = OsString::from(path);
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().map(PathBuf::from).expect("missing src");
    let name = args.next().expect("missing name");

    let reader = File::open(&src)
        .map(BufReader::new)
        .map(fasta::Reader::new)?;

    let fasta_src = push_ext(src, "fai");
    let index = fai::read(fasta_src)?;

    let adapter = IndexedReader::new(reader, index);
    let repository = fasta::Repository::new(adapter);

    let sequence = repository
        .get(&name)
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid name"))?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = fasta::Writer::new(handle);

    let record = fasta::Record::new(Definition::new(name, None), sequence);
    writer.write_record(&record)?;

    Ok(())
}
