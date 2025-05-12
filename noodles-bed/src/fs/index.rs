use std::{fs::File, io, path::Path};

use bstr::ByteSlice;
use noodles_bgzf as bgzf;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_tabix as tabix;

use crate::{io::Reader, Record};

/// Indexes a bgzipped-compressed BED file.
///
/// The input must be coordinate-sorted.
///
/// See also [`tabix::fs::write`] to write the resulting [`tabix::Index`] to a file.
///
/// # Examples
///
/// ```no_run
/// use noodles_bed as bed;
/// let _index = bed::fs::index("sample.bed.gz")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn index<P>(src: P) -> io::Result<tabix::Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(bgzf::io::Reader::new)
        .map(Reader::new)?;

    index_inner(&mut reader)
}

fn index_inner<R>(reader: &mut Reader<3, R>) -> io::Result<tabix::Index>
where
    R: bgzf::io::BufRead,
{
    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::binning_index::index::header::Builder::bed().build());

    let mut record = Record::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let reference_sequence_name = record
            .reference_sequence_name()
            .to_str()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let start = record.feature_start()?;

        let end = record
            .feature_end()
            .transpose()?
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing feature end"))?;

        indexer.add_record(reference_sequence_name, start, end, chunk)?;

        start_position = end_position;
    }

    Ok(indexer.build())
}
