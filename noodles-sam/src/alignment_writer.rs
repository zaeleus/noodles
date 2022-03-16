use std::io;

use super::{AlignmentRecord, Header};

/// An alignment format writer.
pub trait AlignmentWriter {
    /// Writes a SAM header.
    fn write_alignment_header(&mut self, header: &Header) -> io::Result<()>;

    /// Writes an alignment record.
    fn write_alignment_record(
        &mut self,
        header: &Header,
        record: &dyn AlignmentRecord,
    ) -> io::Result<()>;
}
