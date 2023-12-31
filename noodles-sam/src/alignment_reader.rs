use std::io;

use super::{alignment::RecordBuf, Header};

/// An alignment format reader.
pub trait AlignmentReader<R> {
    /// Reads a SAM header.
    fn read_alignment_header(&mut self) -> io::Result<Header>;

    /// Returns an iterator over records.
    fn alignment_records<'a>(
        &'a mut self,
        header: &'a Header,
    ) -> Box<dyn Iterator<Item = io::Result<RecordBuf>> + 'a>;
}
