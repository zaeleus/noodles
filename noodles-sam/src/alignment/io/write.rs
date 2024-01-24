use std::io;

use crate::{alignment::Record, Header};

/// An alignment writer.
///
/// A call to [`Self::finish`] must be made before the writer is dropped.
pub trait Write {
    /// Writes a SAM header.
    fn write_alignment_header(&mut self, header: &Header) -> io::Result<()>;

    /// Writes an alignment record.
    fn write_alignment_record(&mut self, header: &Header, record: &dyn Record) -> io::Result<()>;

    /// Shuts down an alignment writer.
    fn finish(&mut self, header: &Header) -> io::Result<()>;
}
