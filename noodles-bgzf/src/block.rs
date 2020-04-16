use std::io::Cursor;

use super::VirtualPosition;

#[derive(Debug, Default)]
pub struct Block {
    position: u64,
    data: Cursor<Vec<u8>>,
}

impl Block {
    pub fn position(&self) -> u64 {
        self.position
    }

    pub fn set_position(&mut self, pos: u64) {
        self.position = pos;
    }

    pub fn virtual_position(&self) -> VirtualPosition {
        VirtualPosition::from((self.position, self.data.position()))
    }

    pub fn data(&self) -> &Cursor<Vec<u8>> {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut Cursor<Vec<u8>> {
        &mut self.data
    }
}
