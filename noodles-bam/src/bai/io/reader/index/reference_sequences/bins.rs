use std::io::{self, Read};

use indexmap::IndexMap;
use noodles_csi::binning_index::index::reference_sequence::{Bin, Metadata};

use crate::io::reader::num::read_u32_le;

pub(super) fn read_bins<R>(reader: &mut R) -> io::Result<(IndexMap<usize, Bin>, Option<Metadata>)>
where
    R: Read,
{
    use noodles_csi::io::reader::index::reference_sequences::{bins::read_chunks, read_metadata};

    use crate::bai::DEPTH;

    const METADATA_ID: usize = Bin::metadata_id(DEPTH);

    fn duplicate_bin_error(id: usize) -> io::Result<(IndexMap<usize, Bin>, Option<Metadata>)> {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("duplicate bin ID: {id}"),
        ))
    }

    let n_bin = read_u32_le(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bins = IndexMap::with_capacity(n_bin);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = read_u32_le(reader).and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        if id == METADATA_ID {
            let m =
                read_metadata(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            if metadata.replace(m).is_some() {
                return duplicate_bin_error(id);
            }
        } else {
            let chunks =
                read_chunks(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let bin = Bin::new(chunks);

            if bins.insert(id, bin).is_some() {
                return duplicate_bin_error(id);
            }
        }
    }

    Ok((bins, metadata))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_bins() -> io::Result<()> {
        let data = [
            0x00, 0x00, 0x00, 0x00, // n_bin = 0
        ];
        let mut reader = &data[..];
        let (actual_bins, actual_metadata) = read_bins(&mut reader)?;
        assert!(actual_bins.is_empty());
        assert!(actual_metadata.is_none());

        let data = [
            0x02, 0x00, 0x00, 0x00, // n_bin = 2
            0x00, 0x00, 0x00, 0x00, // bins[0].id = 0
            0x00, 0x00, 0x00, 0x00, // bins[0].n_chunk = 0
            0x4a, 0x92, 0x00, 0x00, // bins[1].id = 37450
            0x02, 0x00, 0x00, 0x00, // bins[1].n_chunk = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.n_unmapped = 0
        ];
        let mut reader = &data[..];
        let (actual_bins, actual_metadata) = read_bins(&mut reader)?;
        assert_eq!(actual_bins.len(), 1);
        assert!(actual_bins.get(&0).is_some());
        assert!(actual_metadata.is_some());

        let data = [
            0x01, 0x00, 0x00, 0x00, // n_bin = 1
        ];
        let mut reader = &data[..];
        assert!(matches!(
            read_bins(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = [
            0x02, 0x00, 0x00, 0x00, // n_bin = 2
            0x00, 0x00, 0x00, 0x00, // bins[0].id = 0
            0x00, 0x00, 0x00, 0x00, // bins[0].n_chunk = 0
            0x00, 0x00, 0x00, 0x00, // bins[1].id = 0
            0x00, 0x00, 0x00, 0x00, // bins[1].n_chunk = 0
        ];
        let mut reader = &data[..];
        assert!(matches!(
            read_bins(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        let data = [
            0x02, 0x00, 0x00, 0x00, // n_bin = 2
            0x4a, 0x92, 0x00, 0x00, // bins[0].id = 37450
            0x02, 0x00, 0x00, 0x00, // bins[0].n_chunk = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[0].metadata.ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[0].metadata.ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[0].metadata.n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[0].metadata.n_unmapped = 0
            0x4a, 0x92, 0x00, 0x00, // bins[1].id = 37450
            0x02, 0x00, 0x00, 0x00, // bins[1].n_chunk = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.n_unmapped = 0
        ];
        let mut reader = &data[..];
        assert!(matches!(
            read_bins(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
