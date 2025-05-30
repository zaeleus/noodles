use std::{fs::File, io::BufReader, path::Path};

use crate::{fai, io::Indexer};

/// Indexes a FASTQ file.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_fastq as fastq;
/// let index = fastq::fs::index("sample.fastq")?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn index<P>(src: P) -> std::io::Result<fai::Index>
where
    P: AsRef<Path>,
{
    let mut indexer = File::open(src).map(BufReader::new).map(Indexer::new)?;
    let mut index = Vec::new();

    while let Some(record) = indexer.index_record()? {
        index.push(record);
    }

    Ok(index)
}
