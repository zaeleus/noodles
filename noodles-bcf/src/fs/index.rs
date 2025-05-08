use std::{fs::File, io, path::Path};

use noodles_bgzf as bgzf;
use noodles_csi::{
    self as csi,
    binning_index::{index::reference_sequence::bin::Chunk, Indexer},
};
use noodles_vcf::variant::Record as _;

use crate::{io::Reader, Record};

/// Indexes a BCF file.
///
/// The input must be coordinate-sorted.
///
/// See also [`csi::fs::write`] to write the resulting [`csi::Index`].
///
/// # Examples
///
/// ```no_run
/// use noodles_bcf as bcf;
/// let _index = bcf::fs::index("sample.bcf")?;
/// Ok::<_, std::io::Error>(())
/// ```
pub fn index<P>(src: P) -> io::Result<csi::Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    index_inner(&mut reader)
}

fn index_inner<R>(reader: &mut Reader<R>) -> io::Result<csi::Index>
where
    R: bgzf::io::Read,
{
    let header = reader.read_header()?;
    let mut indexer = Indexer::default();

    let mut record = Record::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let reference_sequence_id = record.reference_sequence_id()?;

        let start = record
            .variant_start()
            .transpose()?
            .expect("missing variant start");

        let end = record.variant_end(&header)?;

        indexer.add_record(Some((reference_sequence_id, start, end, true)), chunk)?;

        start_position = end_position;
    }

    Ok(indexer.build(header.contigs().len()))
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;
    use noodles_csi::{binning_index::ReferenceSequence as _, BinningIndex};
    use noodles_vcf::{
        self as vcf,
        variant::{io::Write, RecordBuf},
    };

    use super::*;

    #[test]
    fn test_index() -> io::Result<()> {
        let mut writer = crate::io::Writer::new(Vec::new());

        let header = vcf::Header::builder()
            .add_contig("sq0", Default::default())
            .build();

        writer.write_header(&header)?;

        let record = RecordBuf::builder()
            .set_reference_sequence_name("sq0")
            .set_variant_start(Position::MIN)
            .set_reference_bases("N")
            .build();

        writer.write_variant_record(&header, &record)?;
        writer.try_finish()?;

        let src = writer.into_inner().into_inner();
        let mut reader = Reader::new(&src[..]);

        let index = index_inner(&mut reader)?;

        assert_eq!(index.min_shift(), 14);
        assert_eq!(index.depth(), 5);
        assert!(index.header().is_none());

        let reference_sequences = index.reference_sequences();
        assert_eq!(reference_sequences.len(), 1);

        let sq0 = &reference_sequences[0];
        let bins = sq0.bins();
        assert_eq!(bins.len(), 1);
        assert!(bins.contains_key(&4681));

        assert_eq!(
            sq0.metadata().map(|metadata| (
                metadata.mapped_record_count(),
                metadata.unmapped_record_count()
            )),
            Some((1, 0))
        );

        assert_eq!(index.unplaced_unmapped_record_count(), Some(0));

        Ok(())
    }
}
