pub(crate) mod header;
mod reference_sequences;

use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use self::{header::read_aux, reference_sequences::read_reference_sequences};
use crate::Index;

pub(super) fn read_index<R>(reader: &mut R) -> io::Result<Index>
where
    R: Read,
{
    read_magic(reader)?;

    let min_shift = reader
        .read_i32::<LittleEndian>()
        .and_then(|n| u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))?;

    let depth = reader
        .read_i32::<LittleEndian>()
        .and_then(|n| u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))?;

    let header = read_aux(reader)?;

    let reference_sequences = read_reference_sequences(reader, depth)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let n_no_coor = read_unplaced_unmapped_record_count(reader)?;

    let mut builder = Index::builder()
        .set_min_shift(min_shift)
        .set_depth(depth)
        .set_reference_sequences(reference_sequences);

    if let Some(hdr) = header {
        builder = builder.set_header(hdr);
    }

    if let Some(n_no_coor) = n_no_coor {
        builder = builder.set_unplaced_unmapped_record_count(n_no_coor);
    }

    Ok(builder.build())
}

fn read_magic<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    use crate::io::MAGIC_NUMBER;

    let mut magic = [0; 4];
    reader.read_exact(&mut magic)?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid CSI file format",
        ))
    }
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
    fn test_read_magic() {
        let data = b"CSI\x01";
        let mut reader = &data[..];
        assert!(read_magic(&mut reader).is_ok());
    }

    #[test]
    fn test_read_magic_with_invalid_magic_number() {
        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"CSI";
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
