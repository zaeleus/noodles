use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::codecs::rans_4x8::Order;

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<(Order, usize, usize)>
where
    R: Read,
{
    let order = reader.read_u8().and_then(|order| {
        Order::try_from(order).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let compressed_len = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let data_len = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    Ok((order, compressed_len, data_len))
}
