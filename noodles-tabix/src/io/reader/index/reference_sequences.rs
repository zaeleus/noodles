mod bins;
mod intervals;

use std::io::{self, Read};

use noodles_csi::binning_index::index::{
    ReferenceSequence, reference_sequence::index::LinearIndex,
};

use self::{bins::read_bins, intervals::read_intervals};

pub(super) fn read_reference_sequences<R>(
    reader: &mut R,
    reference_sequence_count: usize,
) -> io::Result<Vec<ReferenceSequence<LinearIndex>>>
where
    R: Read,
{
    (0..reference_sequence_count)
        .map(|_| read_reference_sequence(reader))
        .collect()
}

fn read_reference_sequence<R>(reader: &mut R) -> io::Result<ReferenceSequence<LinearIndex>>
where
    R: Read,
{
    let (bins, metadata) = read_bins(reader)?;
    let intervals = read_intervals(reader)?;
    Ok(ReferenceSequence::new(bins, intervals, metadata))
}
