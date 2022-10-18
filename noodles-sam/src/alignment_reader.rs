use std::io;

use noodles_fasta as fasta;

use super::{alignment::Record, Header};

/// An alignment format reader.
pub trait AlignmentReader<R> {
    /// Reads a SAM header.
    fn read_alignment_header(&mut self) -> io::Result<Header>;

    /// Returns an iterator over records.
    fn alignment_records<'a>(
        &'a mut self,
        reference_sequence_repository: &'a fasta::Repository,
        header: &'a Header,
    ) -> Box<dyn Iterator<Item = io::Result<Record>> + 'a>;
}
