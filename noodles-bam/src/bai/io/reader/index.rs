mod magic_number;
mod reference_sequences;

use std::io::{self, Read};

use self::{magic_number::read_magic_number, reference_sequences::read_reference_sequences};
use crate::{bai::Index, io::reader::num::read_u64_le};

pub(super) fn read_index<R>(reader: &mut R) -> io::Result<Index>
where
    R: Read,
{
    read_magic_number(reader)?;

    let reference_sequences = read_reference_sequences(reader)?;
    let unplaced_unmapped_record_count = read_unplaced_unmapped_record_count(reader)?;

    let mut builder = Index::builder().set_reference_sequences(reference_sequences);

    if let Some(n) = unplaced_unmapped_record_count {
        builder = builder.set_unplaced_unmapped_record_count(n);
    }

    Ok(builder.build())
}

fn read_unplaced_unmapped_record_count<R>(reader: &mut R) -> io::Result<Option<u64>>
where
    R: Read,
{
    match read_u64_le(reader) {
        Ok(n) => Ok(Some(n)),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
