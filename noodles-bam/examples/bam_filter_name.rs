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
use noodles_sam::record::Name;

fn read_names<P>(src: P) -> io::Result<HashSet<Name>>
where
    P: AsRef<Path>,
{
    let reader = File::open(src).map(BufReader::new)?;
    let mut names = HashSet::new();

    for result in reader.lines() {
        let name = result.and_then(|s| {
            Name::try_from(s.into_bytes())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        names.insert(name);
    }

    Ok(names)
}

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);
    let names_src = args.next().expect("missing names_src");
    let src = args.next().expect("missing src");

    let names = read_names(names_src)?;

    let mut reader = bam::reader::Builder.build_from_path(src)?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = bam::Writer::new(stdout);

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;

        if let Some(name) = record.name() {
            if names.contains(name) {
                writer.write_record(&header, &record)?;
            }
        }
    }

    Ok(())
}
