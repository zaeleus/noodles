use std::{
    io::Cursor,
    ops::{Deref, DerefMut},
};

#[derive(Debug)]
pub struct Block {
    c_offset: u64,
    inner: Cursor<Vec<u8>>,
}

impl Block {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn position(&self) -> u64 {
        self.inner.position()
    }

    pub fn virtual_position(&self) -> u64 {
        self.c_offset() << 16 | self.u_offset()
    }

    pub fn c_offset(&self) -> u64 {
        self.c_offset
    }

    pub fn set_c_offset(&mut self, c_offset: u64) {
        self.c_offset = c_offset;
    }

    pub fn u_offset(&self) -> u64 {
        self.position()
    }
}

impl Default for Block {
    fn default() -> Self {
        Self {
            c_offset: 0,
            inner: Cursor::new(Vec::new()),
        }
    }
}

impl Deref for Block {
    type Target = Cursor<Vec<u8>>;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl DerefMut for Block {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.inner
    }
}
