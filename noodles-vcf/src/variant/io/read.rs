use std::io;

use crate::{variant::Record, Header};

/// A variant format reader.
pub trait Read<R> {
    /// Reads a VCF header.
    fn read_variant_header(&mut self) -> io::Result<Header>;

    /// Returns an iterator over records.
    fn variant_records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<Box<dyn Record>>> + 'r>;
}
