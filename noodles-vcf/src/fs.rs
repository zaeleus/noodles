//! VCF filesystem operations.

use std::{fs::File, io, path::Path};

use noodles_bgzf as bgzf;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_tabix as tabix;

use super::{io::Reader, variant::Record as _, Record};

/// Indexes a bgzipped-compressed VCF file.
///
/// # Examples
///
/// ```no_run
/// use noodles_vcf as vcf;
/// let _index = vcf::fs::index("sample.vcf.gz")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn index<P>(src: P) -> io::Result<tabix::Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bgzf::Reader::new).map(Reader::new)?;
    let header = reader.read_header()?;

    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::binning_index::index::header::Builder::vcf().build());

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
