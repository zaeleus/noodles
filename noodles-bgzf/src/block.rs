use std::{convert::TryFrom, io::Cursor};

use super::{virtual_position, VirtualPosition};

/// A BGZF block.
///
/// A BGZF block is a gzip stream less than 64 KiB and contains an extra field describing the size
/// of the block itself.
#[derive(Debug, Default)]
pub(crate) struct Block {
    position: u64,
    data: Cursor<Vec<u8>>,
}

impl Block {
    /// Sets the position of this block in the compressed stream.
    pub fn set_position(&mut self, pos: u64) {
        self.position = pos;
    }

    /// Returns the virtual position at the current position in the uncompressed data stream.
    pub fn virtual_position(&self) -> VirtualPosition {
        assert!(self.position <= virtual_position::MAX_COMPRESSED_POSITION);
        assert!(self.data.position() <= u64::from(virtual_position::MAX_UNCOMPRESSED_POSITION));

        let uncompressed_position = self.data.position() as u16;

        VirtualPosition::try_from((self.position, uncompressed_position)).unwrap()
    }

    /// Returns a mutable reference to the uncompressed data of this block.
    pub fn data_mut(&mut self) -> &mut Cursor<Vec<u8>> {
        &mut self.data
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_set_position() {
        let mut block = Block::default();
        block.set_position(13);
        assert_eq!(block.position, 13);
    }

    #[test]
    fn test_virtual_position() {
        let mut block = Block::default();
        assert_eq!(block.virtual_position(), VirtualPosition::from(0));

        block.set_position(13);
        assert_eq!(block.virtual_position(), VirtualPosition::from(851968));

        block.data_mut().set_position(2);
        assert_eq!(block.virtual_position(), VirtualPosition::from(851970));
    }

    #[test]
    fn test_data_mut() {
        let mut block = Block::default();
        block.data_mut().get_mut().extend(&[0, 0]);
        assert_eq!(block.data.get_ref().len(), 2);
    }
}
