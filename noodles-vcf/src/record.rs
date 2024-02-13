//! Lazily-evaluated VCF record and fields.

pub(crate) mod fields;
mod filters;
mod ids;
mod info;
mod samples;

use std::fmt;

use self::fields::Fields;
pub use self::{filters::Filters, ids::Ids, info::Info, samples::Samples};

/// An immutable, lazily-evaluated VCF record.
#[derive(Clone, Default, Eq, PartialEq)]
pub struct Record(Fields);

impl Record {
    #[cfg(test)]
    pub(crate) fn fields(&self) -> &Fields {
        &self.0
    }

    pub(crate) fn fields_mut(&mut self) -> &mut Fields {
        &mut self.0
    }

    /// Returns the chromosome.
    pub fn chromosome(&self) -> &str {
        self.0.reference_sequence_name()
    }

    /// Returns the position.
    pub fn position(&self) -> &str {
        self.0.position()
    }

    /// Returns the IDs.
    pub fn ids(&self) -> Ids<'_> {
        self.0.ids()
    }

    /// Returns the reference bases.
    pub fn reference_bases(&self) -> &str {
        self.0.reference_bases()
    }

    /// Returns the alternate bases.
    pub fn alternate_bases(&self) -> &str {
        self.0.alternate_bases()
    }

    /// Returns the quality score.
    pub fn quality_score(&self) -> Option<&str> {
        self.0.quality_score()
    }

    /// Returns the filters.
    pub fn filters(&self) -> Option<Filters<'_>> {
        self.0.filters()
    }

    /// Returns the info.
    pub fn info(&self) -> Info<'_> {
        self.0.info()
    }

    /// Returns the genotypes.
    pub fn genotypes(&self) -> Samples<'_> {
        self.0.samples()
    }
}

impl fmt::Debug for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Record")
            .field("chromosome", &self.chromosome())
            .field("position", &self.position())
            .field("ids", &self.ids())
            .field("reference_bases", &self.reference_bases())
            .field("alternate_bases", &self.alternate_bases())
            .field("quality_score", &self.quality_score())
            .field("filters", &self.filters())
            .field("info", &self.info())
            .field("genotypes", &self.genotypes())
            .finish()
    }
}
