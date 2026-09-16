use std::{fs::File, io, path::Path};

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;

use crate::{Index, Record, index::Indexer, io::Reader, variant::Record as _};

/// Indexes a bgzipped-compressed VCF file.
///
/// This typically returns the index as tabix (TBI); however, if the length of a reference sequence
/// is longer than 2<sup>29</sup>, this returns a coordinate-sorted index (CSI) instead.
///
/// # Examples
///
/// ```no_run
/// use noodles_vcf as vcf;
/// let _index = vcf::fs::index("sample.vcf.gz")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn index<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(bgzf::io::Reader::new)
        .map(Reader::new)?;

    index_inner(&mut reader)
}

fn index_inner<R>(reader: &mut Reader<R>) -> io::Result<Index>
where
    R: bgzf::io::BufRead,
{
    let header = reader.read_header()?;

    let mut builder = Indexer::builder();

    if let Some(max_reference_sequence_length) = header
        .contigs()
        .values()
        .filter_map(|contig| contig.length())
        .max()
    {
        let max_position = Position::try_from(max_reference_sequence_length)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        builder = builder.set_max_position_hint(max_position);
    }

    let mut indexer = builder
        .build()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let mut record = Record::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let reference_sequence_name = record.reference_sequence_name();

        let start = record
            .variant_start()
            .transpose()?
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing position"))?;

        let end = record.variant_end(&header)?;

        indexer.add_record(reference_sequence_name, start, end, chunk)?;

        start_position = end_position;
    }

    Ok(indexer.build())
}
