use std::{
    io::{self, Cursor, Read},
    ops::{Deref, DerefMut},
};

use byteorder::{LittleEndian, ReadBytesExt};

#[derive(Debug)]
pub struct Block(Cursor<Vec<u8>>);

impl Block {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn read_record(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.0.read(buf)
    }

    pub fn is_eof(&self) -> bool {
        self.0.position() >= self.0.get_ref().len() as u64
    }

    pub fn read_block_size(&mut self) -> io::Result<i32> {
        self.0.read_i32::<LittleEndian>()
    }
}

impl Default for Block {
    fn default() -> Self {
        Self(Cursor::new(Vec::new()))
    }
}

impl Deref for Block {
    type Target = Cursor<Vec<u8>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Block {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
