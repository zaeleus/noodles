use std::{
    io::{self, Read},
    num::NonZeroUsize,
};

use bstr::BString;
use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam::header::{
    record::value::{map::ReferenceSequence, Map},
    ReferenceSequences,
};

use crate::io::reader::bytes_with_nul_to_bstring;

pub(super) fn read_reference_sequences<R>(reader: &mut R) -> io::Result<ReferenceSequences>
where
    R: Read,
{
    let n_ref = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut reference_sequences = ReferenceSequences::with_capacity(n_ref);

    for _ in 0..n_ref {
        let (name, reference_sequence) = read_reference_sequence(reader)?;
        reference_sequences.insert(name, reference_sequence);
    }

    Ok(reference_sequences)
}

fn read_reference_sequence<R>(reader: &mut R) -> io::Result<(BString, Map<ReferenceSequence>)>
where
    R: Read,
{
    let l_name = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut c_name = vec![0; l_name];
    reader.read_exact(&mut c_name)?;

    let name = bytes_with_nul_to_bstring(&c_name)?;

    let l_ref = reader.read_u32::<LittleEndian>().and_then(|len| {
        usize::try_from(len)
            .and_then(NonZeroUsize::try_from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let reference_sequence = Map::<ReferenceSequence>::new(l_ref);

    Ok((name, reference_sequence))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_reference_sequences() -> Result<(), Box<dyn std::error::Error>> {
        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        let data = [
            0x01, 0x00, 0x00, 0x00, // n_ref = 1
            0x04, 0x00, 0x00, 0x00, // ref[0].l_name = 4
            0x73, 0x71, 0x30, 0x00, // ref[0].name = "sq0\x00"
            0x08, 0x00, 0x00, 0x00, // ref[0].l_ref = 8
        ];

        let mut reader = &data[..];
        let actual = read_reference_sequences(&mut reader)?;

        let expected: ReferenceSequences =
            [(BString::from("sq0"), Map::<ReferenceSequence>::new(SQ0_LN))]
                .into_iter()
                .collect();

        assert_eq!(actual, expected);

        Ok(())
    }
}
