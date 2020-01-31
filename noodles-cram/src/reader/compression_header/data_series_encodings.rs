use std::io::{self, Read};

use crate::num::read_itf8;

pub fn read_data_series_encodings<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    let data_len = read_itf8(reader)?;
    let mut buf = vec![0; data_len as usize];
    reader.read_exact(&mut buf)
}
