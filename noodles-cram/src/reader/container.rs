use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::{
    container::{Header, ReferenceSequenceId},
    num::{read_itf8, read_ltf8, Itf8},
};

pub fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: Read,
{
    let length = reader.read_i32::<LittleEndian>()?;

    let reference_sequence_id = read_itf8(reader).and_then(|n| {
        ReferenceSequenceId::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let starting_position_on_the_reference = read_itf8(reader)?;
    let alignment_span = read_itf8(reader)?;
    let number_of_records = read_itf8(reader)?;
    let record_counter = read_ltf8(reader)?;
    let bases = read_ltf8(reader)?;
    let number_of_blocks = read_itf8(reader)?;
    let landmarks = read_landmarks(reader)?;
    let crc32 = reader.read_u32::<LittleEndian>()?;

    Ok(Header::new(
        length,
        reference_sequence_id,
        starting_position_on_the_reference,
        alignment_span,
        number_of_records,
        record_counter,
        bases,
        number_of_blocks,
        landmarks,
        crc32,
    ))
}

fn read_landmarks<R>(reader: &mut R) -> io::Result<Vec<Itf8>>
where
    R: Read,
{
    let len = read_itf8(reader).map(|l| l as usize)?;
    let mut buf = Vec::with_capacity(len);

    for _ in 0..len {
        let pos = read_itf8(reader)?;
        buf.push(pos);
    }

    Ok(buf)
}
