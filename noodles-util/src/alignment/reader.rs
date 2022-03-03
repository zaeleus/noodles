mod builder;

pub use self::builder::Builder;

use std::io::{self, Read, Seek};

use noodles_fasta as fasta;
use noodles_sam as sam;

/// An alignment reader.
pub struct Reader {
    inner: Box<dyn sam::AlignmentReader>,
    reference_sequence_repository: fasta::Repository,
}

impl Reader {
    /// Creates an alignment reader builder.
    pub fn builder<R>(inner: R) -> Builder<R>
    where
        R: Read + Seek + 'static,
    {
        Builder::new(inner)
    }

    /// Reads and parses an alignment header.
    pub fn read_header(&mut self) -> io::Result<sam::Header> {
        self.inner.read_alignment_header()
    }

    /// Returns an iterator over records starting from the current stream position.
    pub fn records<'a>(
        &'a mut self,
        header: &'a sam::Header,
    ) -> Box<dyn Iterator<Item = io::Result<sam::Record>> + 'a> {
        self.inner
            .alignment_records(&self.reference_sequence_repository, header)
    }
}
