use std::{
    io::{self, Read},
    num::NonZeroUsize,
};

use bstr::BString;
use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam::header::record::value::{map::ReferenceSequence, Map};

use crate::io::reader::bytes_with_nul_to_bstring;

pub(super) fn read_reference_sequence<R>(
    reader: &mut R,
) -> io::Result<(BString, Map<ReferenceSequence>)>
where
    R: Read,
{
    let name = read_name(reader)?;

    let len = read_length(reader)?;
    let reference_sequence = Map::<ReferenceSequence>::new(len);

    Ok((name, reference_sequence))
}

fn read_name<R>(reader: &mut R) -> io::Result<BString>
where
    R: Read,
{
    let l_name = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut c_name = vec![0; l_name];
    reader.read_exact(&mut c_name)?;

    bytes_with_nul_to_bstring(&c_name)
}

fn read_length<R>(reader: &mut R) -> io::Result<NonZeroUsize>
where
    R: Read,
{
    reader.read_u32::<LittleEndian>().and_then(|len| {
        usize::try_from(len)
            .and_then(NonZeroUsize::try_from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_reference_sequence() -> Result<(), Box<dyn std::error::Error>> {
        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        let data = [
            0x04, 0x00, 0x00, 0x00, // l_name = 4
            0x73, 0x71, 0x30, 0x00, // name = "sq0\x00"
            0x08, 0x00, 0x00, 0x00, // l_ref = 8
        ];

        let mut reader = &data[..];
        let actual = read_reference_sequence(&mut reader)?;
        let expected = (BString::from("sq0"), Map::<ReferenceSequence>::new(SQ0_LN));
        assert_eq!(actual, expected);

        Ok(())
    }
}
