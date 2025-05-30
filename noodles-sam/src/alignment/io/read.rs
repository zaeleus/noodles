use std::io;

use crate::{Header, alignment::Record};

/// An alignment reader.
pub trait Read<R> {
    /// Reads a SAM header.
    fn read_alignment_header(&mut self) -> io::Result<Header>;

    /// Returns an iterator over records.
    fn alignment_records<'a>(
        &'a mut self,
        header: &'a Header,
    ) -> Box<dyn Iterator<Item = io::Result<Box<dyn Record>>> + 'a>;
}
