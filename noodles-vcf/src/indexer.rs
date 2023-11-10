use std::{fs::File, io, path::Path};

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::{self as csi, index::reference_sequence::bin::Chunk};
use noodles_tabix as tabix;

use super::{Reader, Record};

/// Indexes a bgzipped-compressed VCF file.
///
/// ```no_run
/// use noodles_vcf as vcf;
/// let index = vcf::index("sample.vcf.gz")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn index<P>(src: P) -> io::Result<csi::Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bgzf::Reader::new).map(Reader::new)?;
    let header = reader.read_header()?;

    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::index::header::Builder::vcf().build());

    let mut record = Record::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&header, &mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let reference_sequence_name = record.chromosome().to_string();
        let start = Position::try_from(usize::from(record.position()))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let end = record
            .end()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|position| {
                Position::try_from(usize::from(position))
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })?;

        indexer.add_record(&reference_sequence_name, start, end, chunk)?;

        start_position = end_position;
    }

    Ok(indexer.build())
}
