use std::io::{self, Read};

use super::{virtual_position, VirtualPosition};

pub(crate) const MAX_UNCOMPRESSED_DATA_LENGTH: usize = 1 << 16; // bytes

/// A BGZF block.
///
/// A BGZF block is a gzip stream less than 64 KiB and contains an extra field describing the size
/// of the block itself.
#[derive(Debug, Default)]
pub struct Block {
    data: Vec<u8>,
    cpos: u64,
    clen: u64,
    upos: u32,
}

impl Block {
    /// Returns a mutable reference to the uncompressed data of this block.
    pub fn data_mut(&mut self) -> &mut Vec<u8> {
        &mut self.data
    }

    /// Returns the unconsumed part of the current block.
    pub fn fill_buf(&self) -> &[u8] {
        let start = self.upos as usize;
        &self.data[start..]
    }

    /// Returns whether the cursor is at the end of the uncompressed data.
    pub fn is_eof(&self) -> bool {
        self.upos >= self.ulen()
    }

    /// Sets the compressed data length.
    pub fn set_clen(&mut self, clen: u64) {
        self.clen = clen;
    }

    /// Returns the compressed data length.
    #[cfg(feature = "async")]
    pub fn clen(&self) -> u64 {
        self.clen
    }

    /// Sets the position of this block in the compressed stream.
    pub fn set_cpos(&mut self, cpos: u64) {
        self.cpos = cpos;
    }

    /// Sets the position of this block in the uncompressed stream.
    pub fn set_upos(&mut self, upos: u32) {
        self.upos = upos;
    }

    /// Returns the virtual position at the current position in the uncompressed data stream.
    pub fn virtual_position(&self) -> VirtualPosition {
        if self.is_eof() {
            let next_cpos = self.cpos + self.clen;
            assert!(next_cpos <= virtual_position::MAX_COMPRESSED_POSITION);
            VirtualPosition::try_from((next_cpos, 0)).unwrap()
        } else {
            assert!(self.cpos <= virtual_position::MAX_COMPRESSED_POSITION);
            assert!(self.upos <= u32::from(virtual_position::MAX_UNCOMPRESSED_POSITION));
            VirtualPosition::try_from((self.cpos, self.upos as u16)).unwrap()
        }
    }

    /// Returns the uncompressed data length.
    pub fn ulen(&self) -> u32 {
        self.data.len() as u32
    }

    /// Returns the position in the uncompressed data.
    pub fn upos(&self) -> u32 {
        self.upos
    }
}

impl Read for Block {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let upos = self.upos as usize;
        let read_len = (&self.data[upos..]).read(buf)?;
        self.upos += read_len as u32;
        Ok(read_len)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_set_cpos() {
        let mut block = Block::default();
        block.set_cpos(13);
        assert_eq!(block.cpos, 13);
    }

    #[test]
    fn test_set_clen() {
        let mut block = Block::default();
        block.set_clen(8);
        assert_eq!(block.clen, 8);
    }

    #[test]
    fn test_virtual_position() {
        let mut block = Block {
            data: vec![0; 4],
            cpos: 0,
            clen: 8,
            upos: 0,
        };

        assert_eq!(block.virtual_position(), VirtualPosition::from(0));

        block.set_cpos(13);
        assert_eq!(block.virtual_position(), VirtualPosition::from(851968));

        block.set_upos(2);
        assert_eq!(block.virtual_position(), VirtualPosition::from(851970));

        block.set_upos(4);
        assert_eq!(block.virtual_position(), VirtualPosition::from(1376256));
    }

    #[test]
    fn test_data_mut() {
        let mut block = Block::default();
        block.data_mut().extend(&[0, 0]);
        assert_eq!(block.data.len(), 2);
    }

    #[test]
    fn test_is_eof() {
        let mut block = Block {
            data: vec![0; 4],
            cpos: 0,
            clen: 8,
            upos: 0,
        };

        assert!(!block.is_eof());

        block.set_upos(4);

        assert!(block.is_eof());
    }
}
