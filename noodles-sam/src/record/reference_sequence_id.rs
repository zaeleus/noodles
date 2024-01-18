use std::io;

use super::ReferenceSequenceName;
use crate::{alignment::record::fields::ReferenceSequenceId as _, Header};

/// A raw SAM record reference sequence ID.
#[derive(Debug, Eq, PartialEq)]
pub struct ReferenceSequenceId<'h, 'n> {
    header: &'h Header,
    reference_sequence_name: ReferenceSequenceName<'n>,
}

impl<'h, 'n> ReferenceSequenceId<'h, 'n> {
    pub(super) fn new(
        header: &'h Header,
        reference_sequence_name: ReferenceSequenceName<'n>,
    ) -> Self {
        Self {
            header,
            reference_sequence_name,
        }
    }
}

impl<'h, 'n> crate::alignment::record::fields::ReferenceSequenceId for ReferenceSequenceId<'h, 'n> {
    fn try_to_usize(&self) -> io::Result<usize> {
        self.header
            .reference_sequences()
            .get_index_of(self.reference_sequence_name.as_ref())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid reference sequence name",
                )
            })
    }
}

impl<'h, 'n> TryFrom<ReferenceSequenceId<'h, 'n>> for usize {
    type Error = io::Error;

    fn try_from(reference_sequence_id: ReferenceSequenceId<'h, 'n>) -> Result<Self, Self::Error> {
        reference_sequence_id.try_to_usize()
    }
}
