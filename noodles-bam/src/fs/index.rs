use std::{io, path::Path};

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::binning_index::{Indexer, index::reference_sequence::bin::Chunk};
use noodles_sam::{
    self as sam,
    alignment::Record as _,
    header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER},
};

use crate::{Record, bai, io::Reader};

/// Indexes a BAM file.
///
/// The input must be coordinate-sorted and marked as such in the SAM header, i.e.,
/// `SO:coordinate`.
///
/// See also [`bai::fs::write`] to write the resulting [`bai::Index`] to a file.
///
/// # Examples
///
/// ```no_run
/// use noodles_bam as bam;
/// let _index = bam::fs::index("sample.bam")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn index<P>(src: P) -> io::Result<bai::Index>
where
    P: AsRef<Path>,
{
    let mut reader = super::open(src)?;
    index_inner(&mut reader)
}

fn index_inner<R>(reader: &mut Reader<R>) -> io::Result<bai::Index>
where
    R: bgzf::io::Read,
{
    let header = reader.read_header()?;

    if !is_coordinate_sorted(&header) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid sort order: expected {:?}, got {:?}",
                Some(COORDINATE),
                header
                    .header()
                    .and_then(|hdr| hdr.other_fields().get(&SORT_ORDER))
            ),
        ));
    }

    let mut record = Record::default();

    let mut builder = Indexer::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let alignment_context = match alignment_context(&record)? {
            (Some(id), Some(start), Some(end)) => {
                let is_mapped = !record.flags().is_unmapped();
                Some((id, start, end, is_mapped))
            }
            _ => None,
        };

        builder.add_record(alignment_context, chunk)?;

        start_position = end_position;
    }

    Ok(builder.build(header.reference_sequences().len()))
}

fn is_coordinate_sorted(header: &sam::Header) -> bool {
    header
        .header()
        .and_then(|hdr| hdr.other_fields().get(&SORT_ORDER))
        .map(|sort_order| sort_order == COORDINATE)
        .unwrap_or_default()
}

fn alignment_context(
    record: &Record,
) -> io::Result<(Option<usize>, Option<Position>, Option<Position>)> {
    Ok((
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ))
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use bstr::BString;
    use noodles_csi::{BinningIndex, binning_index::ReferenceSequence as _};
    use noodles_sam::{
        alignment::{RecordBuf, io::Write, record::Flags},
        header::record::value::{
            Map,
            map::{self, ReferenceSequence, builder::BuildError},
        },
    };

    use super::*;

    #[test]
    fn test_index() -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = crate::io::Writer::new(Vec::new());

        let header = sam::Header::builder()
            .set_header(
                Map::<map::Header>::builder()
                    .insert(SORT_ORDER, COORDINATE)
                    .build()?,
            )
            .set_reference_sequences(
                [(
                    BString::from("sq0"),
                    Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
                )]
                .into_iter()
                .collect(),
            )
            .build();

        writer.write_header(&header)?;

        let record = RecordBuf::builder()
            .set_flags(Flags::default())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::MIN)
            .build();

        writer.write_alignment_record(&header, &record)?;

        let record = Record::default();
        writer.write_record(&header, &record)?;

        writer.try_finish()?;

        let data = writer.into_inner().into_inner();

        let mut reader = Reader::new(&data[..]);
        let index = index_inner(&mut reader)?;

        assert_eq!(index.min_shift(), 14);
        assert_eq!(index.depth(), 5);
        assert!(index.header().is_none());

        let reference_sequences = index.reference_sequences();
        assert_eq!(reference_sequences.len(), 1);

        let sq0 = &reference_sequences[0];
        let bins = sq0.bins();
        assert_eq!(bins.len(), 1);
        assert!(bins.get(&4681).is_some());

        assert_eq!(
            sq0.metadata().map(|metadata| (
                metadata.mapped_record_count(),
                metadata.unmapped_record_count()
            )),
            Some((1, 0))
        );

        assert_eq!(index.unplaced_unmapped_record_count(), Some(1));

        Ok(())
    }

    #[test]
    fn test_is_coordinate_sorted() -> Result<(), BuildError> {
        let header = sam::Header::default();
        assert!(!is_coordinate_sorted(&header));

        let header = sam::Header::builder()
            .set_header(
                Map::<map::Header>::builder()
                    .insert(SORT_ORDER, COORDINATE)
                    .build()?,
            )
            .build();

        assert!(is_coordinate_sorted(&header));

        Ok(())
    }
}
