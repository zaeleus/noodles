use std::{convert::TryFrom, io::Cursor};

use super::{virtual_position, VirtualPosition};

/// A BGZF block.
///
/// A BGZF block is a gzip stream less than 64 KiB and contains an extra field describing the size
/// of the block itself.
#[derive(Debug, Default)]
pub(crate) struct Block {
    position: u64,
    len: u64,
    data: Cursor<Vec<u8>>,
}

impl Block {
    /// Sets the position of this block in the compressed stream.
    pub fn set_position(&mut self, pos: u64) {
        self.position = pos;
    }

    /// Sets the compressed data length.
    pub fn set_len(&mut self, len: u64) {
        self.len = len;
    }

    /// Returns the virtual position at the current position in the uncompressed data stream.
    pub fn virtual_position(&self) -> VirtualPosition {
        if self.is_eof() {
            let next_compressed_position = self.position + self.len;
            assert!(next_compressed_position <= virtual_position::MAX_COMPRESSED_POSITION);
            VirtualPosition::try_from((next_compressed_position, 0)).unwrap()
        } else {
            assert!(self.position <= virtual_position::MAX_COMPRESSED_POSITION);
            assert!(self.data.position() <= u64::from(virtual_position::MAX_UNCOMPRESSED_POSITION));
            let uncompressed_position = self.data.position() as u16;
            VirtualPosition::try_from((self.position, uncompressed_position)).unwrap()
        }
    }

    /// Returns a mutable reference to the uncompressed data of this block.
    pub fn data_mut(&mut self) -> &mut Cursor<Vec<u8>> {
        &mut self.data
    }

    /// Returns whether the cursor is at the end of the uncompressed data.
    fn is_eof(&self) -> bool {
        let len = self.data.get_ref().len() as u64;
        self.data.position() >= len
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
    fn test_set_len() {
        let mut block = Block::default();
        block.set_len(8);
        assert_eq!(block.len, 8);
    }

    #[test]
    fn test_virtual_position() {
        let mut block = Block {
            position: 0,
            len: 8,
            data: Cursor::new(vec![0; 4]),
        };

        assert_eq!(block.virtual_position(), VirtualPosition::from(0));

        block.set_position(13);
        assert_eq!(block.virtual_position(), VirtualPosition::from(851968));

        block.data_mut().set_position(2);
        assert_eq!(block.virtual_position(), VirtualPosition::from(851970));

        block.data_mut().set_position(4);
        assert_eq!(block.virtual_position(), VirtualPosition::from(1376256));
    }

    #[test]
    fn test_data_mut() {
        let mut block = Block::default();
        block.data_mut().get_mut().extend(&[0, 0]);
        assert_eq!(block.data.get_ref().len(), 2);
    }

    #[test]
    fn test_is_eof() {
        let mut block = Block {
            position: 0,
            len: 8,
            data: Cursor::new(vec![0; 4]),
        };

        assert!(!block.is_eof());

        block.data_mut().set_position(4);

        assert!(block.is_eof());
    }
}
