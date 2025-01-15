use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use noodles_csi::{
    self as csi,
    binning_index::index::{
        reference_sequence::{index::LinearIndex, Bin, Metadata},
        ReferenceSequence,
    },
    io::reader::index::read_header,
};

use crate::{Index, MAGIC_NUMBER};

pub(super) fn read_index<R>(reader: &mut R) -> io::Result<Index>
where
    R: Read,
{
    read_magic(reader)?;

    let n_ref = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let header = read_header(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let references = read_references(reader, n_ref)?;
    let n_no_coor = read_unplaced_unmapped_record_count(reader)?;

    let mut builder = Index::builder()
        .set_header(header)
        .set_reference_sequences(references);

    if let Some(unplaced_unmapped_record_count) = n_no_coor {
        builder = builder.set_unplaced_unmapped_record_count(unplaced_unmapped_record_count);
    }

    Ok(builder.build())
}

fn read_magic<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    let mut magic = [0; MAGIC_NUMBER.len()];
    reader.read_exact(&mut magic)?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid tabix header",
        ))
    }
}

fn read_references<R>(reader: &mut R, len: usize) -> io::Result<Vec<ReferenceSequence<LinearIndex>>>
where
    R: Read,
{
    let mut references = Vec::with_capacity(len);

    for _ in 0..len {
        let (bins, metadata) = read_bins(reader)?;
        let intervals = read_intervals(reader)?;
        references.push(ReferenceSequence::new(bins, intervals, metadata));
    }

    Ok(references)
}

fn read_bins<R>(reader: &mut R) -> io::Result<(IndexMap<usize, Bin>, Option<Metadata>)>
where
    R: Read,
{
    use csi::io::reader::index::reference_sequences::{bins::read_chunks, read_metadata};

    use crate::index::DEPTH;

    const METADATA_ID: usize = Bin::metadata_id(DEPTH);

    let n_bin = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bins = IndexMap::with_capacity(n_bin);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32::<LittleEndian>().and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let is_duplicate = if id == METADATA_ID {
            let m =
                read_metadata(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            metadata.replace(m).is_some()
        } else {
            let chunks =
                read_chunks(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let bin = Bin::new(chunks);

            bins.insert(id, bin).is_some()
        };

        if is_duplicate {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("duplicate bin ID: {id}"),
            ));
        }
    }

    Ok((bins, metadata))
}

fn read_intervals<R>(reader: &mut R) -> io::Result<Vec<bgzf::VirtualPosition>>
where
    R: Read,
{
    let n_intv = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut intervals = Vec::with_capacity(n_intv);

    for _ in 0..n_intv {
        let ioff = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        intervals.push(ioff);
    }

    Ok(intervals)
}

fn read_unplaced_unmapped_record_count<R>(reader: &mut R) -> io::Result<Option<u64>>
where
    R: Read,
{
    match reader.read_u64::<LittleEndian>() {
        Ok(n) => Ok(Some(n)),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_magic_with_invalid_magic_number() {
        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"TBI";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"MThd";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_read_unplaced_unmapped_record_count() -> io::Result<()> {
        let data = [];
        let mut reader = &data[..];
        assert_eq!(read_unplaced_unmapped_record_count(&mut reader)?, None);

        let data = [0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_unplaced_unmapped_record_count(&mut reader)?, Some(8));

        Ok(())
    }
}
