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
use noodles_sam::record::ReadName;

fn read_read_names<P>(src: P) -> io::Result<HashSet<ReadName>>
where
    P: AsRef<Path>,
{
    let reader = File::open(src).map(BufReader::new)?;
    let mut read_names = HashSet::new();

    for result in reader.lines() {
        let read_name = result.and_then(|s| {
            ReadName::try_from(s.into_bytes())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        read_names.insert(read_name);
    }

    Ok(read_names)
}

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);
    let read_names_src = args.next().expect("missing read_names_src");
    let src = args.next().expect("missing src");

    let read_names = read_read_names(read_names_src)?;

    let mut reader = bam::reader::Builder.build_from_path(src)?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = bam::Writer::new(stdout);

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;

        if let Some(read_name) = record.read_name() {
            if read_names.contains(read_name) {
                writer.write_record(&header, &record)?;
            }
        }
    }

    Ok(())
}
