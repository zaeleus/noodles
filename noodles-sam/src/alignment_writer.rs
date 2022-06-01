use std::io;

use super::{alignment::Record, Header};

/// An alignment format writer.
///
/// A call to [`finish`] must be made before the writer is dropped.
pub trait AlignmentWriter {
    /// Writes a SAM header.
    fn write_alignment_header(&mut self, header: &Header) -> io::Result<()>;

    /// Writes an alignment record.
    fn write_alignment_record(&mut self, header: &Header, record: &Record) -> io::Result<()>;

    /// Shuts down an alignment format writer.
    fn finish(&mut self, header: &Header) -> io::Result<()>;
}
