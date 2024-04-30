use std::io;

use crate::VirtualPosition;

/// A BGZF reader.
pub trait Read: io::Read {
    /// Returns the current virtual position.
    fn virtual_position(&self) -> VirtualPosition;
}
