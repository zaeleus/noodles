use std::{io, str};

use super::ReferenceSequenceName;
use crate::Header;

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

impl<'h, 'n> TryFrom<ReferenceSequenceId<'h, 'n>> for usize {
    type Error = io::Error;

    fn try_from(
        ReferenceSequenceId {
            header,
            reference_sequence_name,
        }: ReferenceSequenceId<'h, 'n>,
    ) -> Result<Self, Self::Error> {
        let name = str::from_utf8(reference_sequence_name.as_ref())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        header
            .reference_sequences()
            .get_index_of(name)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid reference sequence name",
                )
            })
    }
}
