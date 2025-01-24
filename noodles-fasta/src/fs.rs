//! FASTA filesystem operations.

use std::{fs::File, io, path::Path};

use super::{fai, io::Indexer};

/// Indexes a FASTA file.
///
/// # Examples
///
/// ```no_run
/// use noodles_fasta as fasta;
/// let index = fasta::fs::index("reference.fa")?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn index<P>(src: P) -> io::Result<fai::Index>
where
    P: AsRef<Path>,
{
    let mut indexer = File::open(src).map(io::BufReader::new).map(Indexer::new)?;
    let mut records = Vec::new();

    while let Some(record) = indexer.index_record()? {
        records.push(record);
    }

    Ok(fai::Index::from(records))
}
