use std::{fs::File, io, num::NonZero, path::Path};

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;

use crate::{Index, Record, header::Contigs, index::Indexer, io::Reader, variant::Record as _};

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

    if let Some(max_len) = max_reference_sequence_length(header.contigs())? {
        // SAFETY: `max_len` is nonzero.
        let max_position = Position::new(max_len.get()).unwrap();
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

        indexer
            .add_record(reference_sequence_name, start, end, chunk)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        start_position = end_position;
    }

    Ok(indexer.build())
}

fn max_reference_sequence_length(contigs: &Contigs) -> io::Result<Option<NonZero<usize>>> {
    let mut max_len = 0;

    for contig in contigs.values() {
        if let Some(len) = contig.length() {
            if len == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid reference sequence length",
                ));
            } else {
                max_len = max_len.max(len);
            }
        }
    }

    Ok(NonZero::new(max_len))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::header::record::value::{Map, map::Contig};

    #[test]
    fn test_max_reference_sequence_length() -> Result<(), Box<dyn std::error::Error>> {
        let contigs = Contigs::new();
        assert!(max_reference_sequence_length(&contigs)?.is_none());

        let contigs = [(
            String::from("sq0"),
            Map::<Contig>::builder().set_length(5).build()?,
        )]
        .into_iter()
        .collect();
        assert_eq!(max_reference_sequence_length(&contigs)?, NonZero::new(5));

        let contigs = [
            (
                String::from("sq0"),
                Map::<Contig>::builder().set_length(5).build()?,
            ),
            (
                String::from("sq1"),
                Map::<Contig>::builder().set_length(13).build()?,
            ),
            (
                String::from("sq2"),
                Map::<Contig>::builder().set_length(8).build()?,
            ),
        ]
        .into_iter()
        .collect();
        assert_eq!(max_reference_sequence_length(&contigs)?, NonZero::new(13));

        let contigs = [(
            String::from("sq0"),
            Map::<Contig>::builder().set_length(0).build()?,
        )]
        .into_iter()
        .collect();
        assert!(matches!(
            max_reference_sequence_length(&contigs),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
