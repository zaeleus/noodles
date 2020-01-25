use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::{
    container::Header,
    num::{read_itf8, read_ltf8, Itf8, Ltf8},
};

pub fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: Read,
{
    let length = read_length(reader)?;
    let reference_sequence_id = read_reference_sequence_id(reader)?;
    let start_pos = read_starting_position_on_the_reference(reader)?;
    let alignment_span = read_alignment_span(reader)?;
    let n_records = read_number_of_records(reader)?;
    let record_counter = read_record_counter(reader)?;
    let bases = read_bases(reader)?;
    let n_blocks = read_number_of_blocks(reader)?;
    let landmarks = read_landmarks(reader)?;
    let _crc = read_crc32(reader)?;

    Ok(Header::new(
        length,
        reference_sequence_id,
        start_pos,
        alignment_span,
        n_records,
        record_counter,
        bases,
        n_blocks,
        landmarks,
    ))
}

fn read_length<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    reader.read_i32::<LittleEndian>()
}

fn read_reference_sequence_id<R>(reader: &mut R) -> io::Result<Itf8>
where
    R: Read,
{
    read_itf8(reader)
}

fn read_starting_position_on_the_reference<R>(reader: &mut R) -> io::Result<Itf8>
where
    R: Read,
{
    read_itf8(reader)
}

fn read_alignment_span<R>(reader: &mut R) -> io::Result<Itf8>
where
    R: Read,
{
    read_itf8(reader)
}

fn read_number_of_records<R>(reader: &mut R) -> io::Result<Itf8>
where
    R: Read,
{
    read_itf8(reader)
}

fn read_record_counter<R>(reader: &mut R) -> io::Result<Ltf8>
where
    R: Read,
{
    read_ltf8(reader)
}

fn read_bases<R>(reader: &mut R) -> io::Result<Ltf8>
where
    R: Read,
{
    read_ltf8(reader)
}

fn read_number_of_blocks<R>(reader: &mut R) -> io::Result<Itf8>
where
    R: Read,
{
    read_itf8(reader)
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

fn read_crc32<R>(reader: &mut R) -> io::Result<[u8; 4]>
where
    R: Read,
{
    let mut buf = [0; 4];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}
