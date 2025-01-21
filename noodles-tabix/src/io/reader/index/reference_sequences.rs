mod bins;
mod intervals;

use std::io::{self, Read};

use noodles_csi::binning_index::index::{
    reference_sequence::index::LinearIndex, ReferenceSequence,
};

use self::{bins::read_bins, intervals::read_intervals};

pub(super) fn read_reference_sequences<R>(
    reader: &mut R,
    reference_sequence_count: usize,
) -> io::Result<Vec<ReferenceSequence<LinearIndex>>>
where
    R: Read,
{
    let mut reference_sequences = Vec::with_capacity(reference_sequence_count);

    for _ in 0..reference_sequence_count {
        let reference_sequence = read_reference_sequence(reader)?;
        reference_sequences.push(reference_sequence);
    }

    Ok(reference_sequences)
}

fn read_reference_sequence<R>(reader: &mut R) -> io::Result<ReferenceSequence<LinearIndex>>
where
    R: Read,
{
    let (bins, metadata) = read_bins(reader)?;
    let intervals = read_intervals(reader)?;
    Ok(ReferenceSequence::new(bins, intervals, metadata))
}
