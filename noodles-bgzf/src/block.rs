use std::io::Cursor;

use super::{virtual_position, VirtualPosition};

/// A BGZF block.
///
/// A BGZF block is a gzip stream less than 64 KiB and contains an extra field describing the size
/// of the block itself.
#[derive(Debug, Default)]
pub struct Block {
    position: u64,
    data: Cursor<Vec<u8>>,
}

impl Block {
    /// Returns the position of this block in the compressed stream.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let block = bgzf::Block::default();
    /// assert_eq!(block.position(), 0);
    /// ```
    pub fn position(&self) -> u64 {
        self.position
    }

    /// Sets the position of this block in the compressed stream.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let mut block = bgzf::Block::default();
    /// block.set_position(13);
    /// assert_eq!(block.position(), 13);
    /// ```
    pub fn set_position(&mut self, pos: u64) {
        self.position = pos;
    }

    /// Returns the virtual position at the current position in the uncompressed data stream.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    ///
    /// let mut block = bgzf::Block::default();
    /// assert_eq!(block.virtual_position(), bgzf::VirtualPosition::from(0));
    ///
    /// block.set_position(13);
    /// assert_eq!(block.virtual_position(), bgzf::VirtualPosition::from(851968));
    ///
    /// block.data_mut().set_position(2);
    /// assert_eq!(block.virtual_position(), bgzf::VirtualPosition::from(851970));
    /// ```
    pub fn virtual_position(&self) -> VirtualPosition {
        assert!(self.data.position() <= u64::from(virtual_position::MAX_UNCOMPRESSED_POSITION));
        let uncompressed_position = self.data.position() as u16;
        VirtualPosition::from((self.position, uncompressed_position))
    }

    /// Returns the uncompressed data of this block.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let block = bgzf::Block::default();
    /// assert!(block.data().get_ref().is_empty());
    /// ```
    pub fn data(&self) -> &Cursor<Vec<u8>> {
        &self.data
    }

    /// Returns a mutable reference to the uncompressed data of this block.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let mut block = bgzf::Block::default();
    /// block.data_mut().get_mut().extend(&[0, 0]);
    /// assert_eq!(block.data().get_ref().len(), 2);
    /// ```
    pub fn data_mut(&mut self) -> &mut Cursor<Vec<u8>> {
        &mut self.data
    }
}
