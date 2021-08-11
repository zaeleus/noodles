//! Filters records in a BAM that match a given list of read names.
//!
//! This is similar to the functionality of
//!
//! ```
//! picard FilterSamReads \
//!     --FILTER includeReadList \
//!     --READ_LIST_FILE <read-names-src> \
//!     --INPUT <src> \
//!     --OUTPUT /dev/stdout
//! ```
//!
//! or
//!
//! `samtools view --qname-file <read-names-src> <src>`.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::{
    collections::HashSet,
    env,
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

use noodles_bam as bam;
use noodles_sam as sam;

fn read_read_names<P>(src: P) -> io::Result<HashSet<String>>
where
    P: AsRef<Path>,
{
    let reader = File::open(src).map(BufReader::new)?;
    let mut read_names = HashSet::new();

    for result in reader.lines() {
        let read_name = result?;
        read_names.insert(read_name);
    }

    Ok(read_names)
}

fn read_read_name(record: &bam::Record) -> io::Result<&str> {
    record
        .read_name()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|c_read_name| {
            c_read_name
                .to_str()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");
    let read_names_src = args.next().expect("missing read_names_src");

    let read_names = read_read_names(read_names_src)?;

    let mut reader = File::open(src).map(bam::Reader::new)?;
    let header: sam::Header = reader.read_header()?.parse()?;
    reader.read_reference_sequences()?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = bam::Writer::new(handle);

    writer.write_header(&header)?;
    writer.write_reference_sequences(header.reference_sequences())?;

    for result in reader.records() {
        let record = result?;

        let read_name = read_read_name(&record)?;

        if read_names.contains(read_name) {
            writer.write_record(&record)?;
        }
    }

    Ok(())
}
