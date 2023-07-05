mod modification;
mod status;
mod strand;
mod unmodified_base;

use self::{
    modification::Modification, status::Status, strand::Strand, unmodified_base::UnmodifiedBase,
};

/// A base modifications group.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Group {
    unmodified_base: UnmodifiedBase,
    strand: Strand,
    modifications: Vec<Modification>,
    status: Option<Status>,
}

impl Group {
    /// Creates a base modifications group.
    fn new(
        unmodified_base: UnmodifiedBase,
        strand: Strand,
        modifications: Vec<Modification>,
        status: Option<Status>,
    ) -> Self {
        Self {
            unmodified_base,
            strand,
            modifications,
            status,
        }
    }

    /// Returns the unmodified base.
    pub fn unmodified_base(&self) -> UnmodifiedBase {
        self.unmodified_base
    }

    /// Returns the strand.
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Returns the modifications.
    pub fn modifications(&self) -> &[Modification] {
        &self.modifications
    }

    /// Returns the status.
    pub fn status(&self) -> Option<Status> {
        self.status
    }
}
