//! Lazily-evaluated VCF record and fields.

mod bounds;
mod filters;
mod genotypes;
mod ids;
mod info;

use std::fmt;

use self::bounds::Bounds;
pub use self::{filters::Filters, genotypes::Genotypes, ids::Ids, info::Info};

const MISSING: &str = ".";

/// An immutable, lazily-evaluated VCF record.
#[derive(Clone, Eq, PartialEq)]
pub struct Record {
    pub(crate) buf: String,
    pub(crate) bounds: Bounds,
}

impl Record {
    /// Returns the chromosome.
    pub fn chromosome(&self) -> &str {
        &self.buf[self.bounds.chromosome_range()]
    }

    /// Returns the position.
    pub fn position(&self) -> &str {
        &self.buf[self.bounds.position_range()]
    }

    /// Returns the IDs.
    pub fn ids(&self) -> Ids<'_> {
        let buf = &self.buf[self.bounds.ids_range()];
        Ids::new(buf)
    }

    /// Returns the reference bases.
    pub fn reference_bases(&self) -> &str {
        &self.buf[self.bounds.reference_bases_range()]
    }

    /// Returns the alternate bases.
    pub fn alternate_bases(&self) -> &str {
        &self.buf[self.bounds.alternate_bases_range()]
    }

    /// Returns the quality score.
    pub fn quality_score(&self) -> Option<&str> {
        match &self.buf[self.bounds.quality_score_range()] {
            MISSING => None,
            buf => Some(buf),
        }
    }

    /// Returns the filters.
    pub fn filters(&self) -> Option<Filters<'_>> {
        let buf = &self.buf[self.bounds.filters_range()];

        match buf {
            MISSING => None,
            _ => Some(Filters::new(buf)),
        }
    }

    /// Returns the info.
    pub fn info(&self) -> Info<'_> {
        let buf = match &self.buf[self.bounds.info_range()] {
            MISSING => "",
            buf => buf,
        };

        Info::new(buf)
    }

    /// Returns the genotypes.
    pub fn genotypes(&self) -> Genotypes<'_> {
        let buf = &self.buf[self.bounds.genotypes_range()];
        Genotypes::new(buf)
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

impl Default for Record {
    fn default() -> Self {
        let buf = String::from("sq01.A....");

        let bounds = Bounds {
            chromosome_end: 3,
            position_end: 4,
            ids_end: 5,
            reference_bases_end: 6,
            alternate_bases_end: 7,
            quality_score_end: 8,
            filters_end: 9,
            info_end: 10,
        };

        Self { buf, bounds }
    }
}
