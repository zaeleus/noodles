mod data;

use self::data::Data;
use super::{virtual_position, VirtualPosition};

/// A BGZF block.
///
/// A BGZF block is a gzip stream less than 64 KiB and contains an extra field describing the size
/// of the block itself.
#[derive(Debug, Default)]
pub struct Block {
    /// The position of the compressed block.
    pos: u64,
    /// The block size (`BSIZE` + 1).
    size: u64,
    /// The uncompressed data (`inflate(CDATA)`).
    data: Data,
}

impl Block {
    pub fn set_position(&mut self, position: u64) {
        self.pos = position;
    }

    pub fn size(&self) -> u64 {
        self.size
    }

    pub fn set_size(&mut self, size: u64) {
        self.size = size;
    }

    /// Returns the virtual position at the current position in the uncompressed data stream.
    pub fn virtual_position(&self) -> VirtualPosition {
        if self.data.has_remaining() {
            assert!(self.pos <= virtual_position::MAX_COMPRESSED_POSITION);
            assert!(
                self.data.position() <= usize::from(virtual_position::MAX_UNCOMPRESSED_POSITION)
            );
            VirtualPosition::try_from((self.pos, self.data.position() as u16)).unwrap()
        } else {
            let next_cpos = self.pos + self.size;
            assert!(next_cpos <= virtual_position::MAX_COMPRESSED_POSITION);
            VirtualPosition::try_from((next_cpos, 0)).unwrap()
        }
    }

    pub fn data(&self) -> &Data {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut Data {
        &mut self.data
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_virtual_position() {
        let mut block = Block::default();
        block.set_size(8);
        block.data_mut().resize(4);

        assert_eq!(block.virtual_position(), VirtualPosition::from(0));

        block.set_position(13);
        assert_eq!(block.virtual_position(), VirtualPosition::from(851968));

        block.data_mut().set_position(2);
        assert_eq!(block.virtual_position(), VirtualPosition::from(851970));

        block.data_mut().set_position(4);
        assert_eq!(block.virtual_position(), VirtualPosition::from(1376256));
    }
}
