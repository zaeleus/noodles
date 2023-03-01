use std::io;

use super::{Header, Record};

/// A variant format reader.
pub trait VariantReader<R> {
    /// Reads a VCF header.
    fn read_variant_header(&mut self) -> io::Result<Header>;

    /// Returns an iterator over records.
    fn variant_records<'a>(
        &'a mut self,
        header: &'a Header,
    ) -> Box<dyn Iterator<Item = io::Result<Record>> + 'a>;
}
