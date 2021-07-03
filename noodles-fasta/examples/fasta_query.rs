//! Queries a FASTA with a given reference sequence name.
//!
//! The input FASTA must have an index in the same directory.
//!
//! The result is similar to the output of `samtools faidx --length 80 <src>
//! <reference-sequence-name>`.

use std::{
    env,
    fs::File,
    io::{self, BufRead, BufReader, Seek, SeekFrom},
    path::PathBuf,
};

use noodles_fasta::{self as fasta, fai};

fn query<R>(
    reader: &mut fasta::Reader<R>,
    index: &[fai::Record],
    reference_sequence_name: &str,
) -> io::Result<fasta::Record>
where
    R: BufRead + Seek,
{
    let index_record = index
        .iter()
        .find(|r| r.reference_sequence_name() == reference_sequence_name)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "invalid reference sequence name: {}",
                    reference_sequence_name
                ),
            )
        })?;

    let pos = index_record.offset();
    reader.seek(SeekFrom::Start(pos))?;

    let definition = fasta::record::Definition::new(reference_sequence_name.into(), None);
    let mut sequence = Vec::new();

    reader.read_sequence(&mut sequence)?;

    Ok(fasta::Record::new(definition, sequence))
}

fn main() -> io::Result<()> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let reference_sequence_name = args.next().expect("missing reference_sequence_name");

    let mut reader = File::open(&src)
        .map(BufReader::new)
        .map(fasta::Reader::new)?;

    let index = fai::read(src.with_extension("fa.fai"))?;

    let record = query(&mut reader, &index, &reference_sequence_name)?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = fasta::Writer::new(handle);

    writer.write_record(&record)?;

    Ok(())
}
