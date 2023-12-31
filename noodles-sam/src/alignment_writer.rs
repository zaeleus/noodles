use std::io;

use super::{alignment::RecordBuf, Header};

/// An alignment format writer.
///
/// A call to [`Self::finish`] must be made before the writer is dropped.
pub trait AlignmentWriter {
    /// Writes a SAM header.
    fn write_alignment_header(&mut self, header: &Header) -> io::Result<()>;

    /// Writes an alignment record.
    fn write_alignment_record(&mut self, header: &Header, record: &RecordBuf) -> io::Result<()>;

    /// Shuts down an alignment format writer.
    fn finish(&mut self, header: &Header) -> io::Result<()>;
}
