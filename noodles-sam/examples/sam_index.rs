//! Builds and writes a SAM index from a SAM file.
//!
//! The input SAM must be bgzip-compressed and coordinate-sorted, i.e., `SO:coordinate`.
//!
//! This writes the output to stdout rather than `<src>.csi`.
//!
//! The output is similar to the output of `samtools index -c <src>`.

use std::{env, fs::File, io};

use noodles_bgzf as bgzf;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_sam::{self as sam, alignment::Record};

fn is_coordinate_sorted(header: &sam::Header) -> bool {
    use sam::header::record::value::map::header::SortOrder;

    header
        .header()
        .and_then(|hdr| hdr.sort_order())
        .map(|sort_order| sort_order == SortOrder::Coordinate)
        .unwrap_or_default()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(bgzf::Reader::new)
        .map(sam::Reader::new)?;

    let header = reader.read_header()?;

    if !is_coordinate_sorted(&header) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "the input SAM must be coordinate-sorted to be indexed",
        )
        .into());
    }

    let mut record = Record::default();

    let mut indexer = csi::binning_index::index::Indexer::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&header, &mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
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

        indexer.add_record(alignment_context, chunk)?;

        start_position = end_position;
    }

    let index = indexer.build(header.reference_sequences().len());

    let stdout = io::stdout().lock();
    let mut writer = csi::Writer::new(stdout);

    writer.write_index(&index)?;

    Ok(())
}
