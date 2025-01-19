mod bins;
mod intervals;

use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_csi::binning_index::index::{
    reference_sequence::index::LinearIndex, ReferenceSequence,
};

use self::{bins::read_bins, intervals::read_intervals};

pub(super) fn read_reference_sequences<R>(
    reader: &mut R,
) -> io::Result<Vec<ReferenceSequence<LinearIndex>>>
where
    R: Read,
{
    let n_ref = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut references = Vec::with_capacity(n_ref);

    for _ in 0..n_ref {
        let (bins, metadata) = read_bins(reader)?;
        let intervals = read_intervals(reader)?;
        references.push(ReferenceSequence::new(bins, intervals, metadata));
    }

    Ok(references)
}
