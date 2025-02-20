use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::codecs::rans_4x8::Order;

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<(Order, usize, usize)>
where
    R: Read,
{
    let order = read_order(reader)?;
    let compressed_size = read_size(reader)?;
    let uncompressed_size = read_size(reader)?;
    Ok((order, compressed_size, uncompressed_size))
}

fn read_order<R>(reader: &mut R) -> io::Result<Order>
where
    R: Read,
{
    reader.read_u8().and_then(|order| {
        Order::try_from(order).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

fn read_size<R>(reader: &mut R) -> io::Result<usize>
where
    R: Read,
{
    reader
        .read_u32::<LittleEndian>()
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}
