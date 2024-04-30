//! Builds and writes a BAM index from a BAM file.
//!
//! The input BAM must be coordinate-sorted, i.e., `SO:coordinate`.
//!
//! This writes the output to stdout rather than `<src>.bai`.
//!
//! The output is similar to the output of `samtools index <src>`.

use std::{env, io};

use noodles_bam::{self as bam, bai};
use noodles_core::Position;
use noodles_csi::binning_index::{index::reference_sequence::bin::Chunk, Indexer};
use noodles_sam::{self as sam, alignment::Record as _};

fn is_coordinate_sorted(header: &sam::Header) -> bool {
    use sam::header::record::value::map::header::{sort_order, tag};

    header
        .header()
        .and_then(|hdr| hdr.other_fields().get(&tag::SORT_ORDER))
        .map(|sort_order| sort_order == sort_order::COORDINATE)
        .unwrap_or_default()
}

fn alignment_context(
    record: &bam::Record,
) -> io::Result<(Option<usize>, Option<Position>, Option<Position>)> {
    Ok((
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ))
}

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bam::io::reader::Builder.build_from_path(src)?;
    let header = reader.read_header()?;

    if !is_coordinate_sorted(&header) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "the input BAM must be coordinate-sorted to be indexed",
        ));
    }

    let mut record = bam::Record::default();

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

    let index = builder.build(header.reference_sequences().len());

    let stdout = io::stdout().lock();
    let mut writer = bai::Writer::new(stdout);
    writer.write_index(&index)?;

    Ok(())
}
