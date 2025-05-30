use std::io::{self, SeekFrom};

use crate::{VirtualPosition, gzi};

/// A seekable BGZF reader.
pub trait Seek {
    /// Seeks the stream to the given virtual position.
    fn seek_to_virtual_position(&mut self, pos: VirtualPosition) -> io::Result<VirtualPosition>;

    /// Seeks the stream to the given position using an index.
    fn seek_with_index(&mut self, index: &gzi::Index, pos: SeekFrom) -> io::Result<u64>;
}
