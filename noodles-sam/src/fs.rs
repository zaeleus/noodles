//! SAM filesystem operations.

use std::{fs::File, io, path::Path};

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::{
    self as csi,
    binning_index::{index::reference_sequence::bin::Chunk, Indexer},
};

use super::{
    alignment::Record as _,
    header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER},
    io::Reader,
    Header, Record,
};

/// Indexes a bgzipped-compressed SAM file.
///
/// The input must be coordinate-sorted and marked as such in the SAM header, i.e.,
/// `SO:coordinate`.
///
/// See also [`csi::fs::write`] to write the resulting [`csi::Index`] to a file.
///
/// # Examples
///
/// ```no_run
/// use noodles_sam as sam;
/// let index = sam::fs::index("sample.sam.gz")?;
/// # Ok::<(), std::io::Error>(())
/// ````
pub fn index<P>(src: P) -> io::Result<csi::Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(bgzf::io::Reader::new)
        .map(Reader::new)?;

    index_inner(&mut reader)
}

fn index_inner<R>(reader: &mut Reader<R>) -> io::Result<csi::Index>
where
    R: bgzf::io::BufRead,
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

    let mut indexer = Indexer::default();

    let mut record = Record::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let alignment_context = match alignment_context(&header, &record)? {
            (Some(id), Some(start), Some(end)) => {
                let is_mapped = !record.flags()?.is_unmapped();
                Some((id, start, end, is_mapped))
            }
            _ => None,
        };

        indexer.add_record(alignment_context, chunk)?;

        start_position = end_position;
    }

    Ok(indexer.build(header.reference_sequences().len()))
}

fn is_coordinate_sorted(header: &Header) -> bool {
    header
        .header()
        .and_then(|hdr| hdr.other_fields().get(&SORT_ORDER))
        .map(|sort_order| sort_order == COORDINATE)
        .unwrap_or_default()
}

fn alignment_context(
    header: &Header,
    record: &Record,
) -> io::Result<(Option<usize>, Option<Position>, Option<Position>)> {
    Ok((
        record.reference_sequence_id(header).transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::header::record::value::{
        map::{self, builder::BuildError},
        Map,
    };

    #[test]
    fn test_is_coordinate_sorted() -> Result<(), BuildError> {
        let header = Header::default();
        assert!(!is_coordinate_sorted(&header));

        let header = Header::builder()
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
