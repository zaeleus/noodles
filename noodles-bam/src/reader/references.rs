use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::Reference;

use super::bytes_with_nul_to_string;

pub struct References<'a, R: Read> {
    reader: &'a mut R,
    i: usize,
    len: usize,
}

impl<'a, R: Read> References<'a, R> {
    pub fn new(reader: &'a mut R, len: usize) -> References<R> {
        References { reader, i: 0, len }
    }
}

impl<'a, R: 'a + Read> Iterator for References<'a, R> {
    type Item = io::Result<Reference>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i >= self.len {
            return None;
        }

        let result = read_reference(self.reader);

        self.i += 1;

        Some(result)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.len, Some(self.len))
    }
}

fn read_reference<R>(reader: &mut R) -> io::Result<Reference>
where
    R: Read,
{
    let l_name = reader.read_i32::<LittleEndian>()?;
    let name = read_name(reader, l_name as usize)?;
    let l_ref = reader.read_i32::<LittleEndian>()?;
    Ok(Reference::new(name, l_ref))
}

fn read_name<R>(reader: &mut R, l_name: usize) -> io::Result<String>
where
    R: Read,
{
    let mut buf = vec![0; l_name];
    reader.read_exact(&mut buf)?;
    bytes_with_nul_to_string(&buf)
}
