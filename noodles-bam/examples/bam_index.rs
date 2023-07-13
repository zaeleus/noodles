//! Builds and writes a BAM index from a BAM file.
//!
//! The input BAM must be coordinate-sorted, i.e., `SO:coordinate`.
//!
//! This writes the output to stdout rather than `<src>.bai`.
//!
//! The output is similar to the output of `samtools index <src>`.

use std::{env, io};

use noodles_bam::{self as bam, bai};
use noodles_csi::{self as csi, index::reference_sequence::bin::Chunk};
use noodles_sam::{self as sam, alignment::Record};

fn is_coordinate_sorted(header: &sam::Header) -> bool {
    use sam::header::record::value::map::header::SortOrder;

    if let Some(hdr) = header.header() {
        if let Some(sort_order) = hdr.sort_order() {
            return sort_order == SortOrder::Coordinate;
        }
    }

    false
}

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bam::reader::Builder.build_from_path(src)?;
    let header = reader.read_header()?;

    if !is_coordinate_sorted(&header) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "the input BAM must be coordinate-sorted to be indexed",
        ));
    }

    let mut record = Record::default();

    let mut builder = csi::index::Indexer::default();
    let mut start_position = reader.virtual_position();

    while reader.read_record(&header, &mut record)? != 0 {
        let end_position = reader.virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let alignment_context = match (
            record.reference_sequence_id(),
            record.alignment_start(),
            record.alignment_end(),
        ) {
            (Some(id), Some(start), Some(end)) => {
                Some((id, start, end, !record.flags().is_unmapped()))
            }
            _ => None,
        };

        builder.add_record(alignment_context, chunk)?;

        start_position = end_position;
    }

    let index = builder.build(header.reference_sequences().len());

    let stdout = io::stdout().lock();
    let mut writer = bai::Writer::new(stdout);

    writer.write_header()?;
    writer.write_index(&index)?;

    Ok(())
}
