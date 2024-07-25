//! VCF record.

mod alternate_bases;
pub(crate) mod fields;
mod filters;
mod ids;
mod info;
pub mod samples;

use std::{fmt, io};

use noodles_core::Position;

use self::fields::Fields;
pub use self::{
    alternate_bases::AlternateBases, filters::Filters, ids::Ids, info::Info, samples::Samples,
};
use super::Header;

/// A VCF record.
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

    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &str {
        self.0.reference_sequence_name()
    }

    /// Returns the variant start position.
    pub fn variant_start(&self) -> Option<io::Result<Position>> {
        self.0.variant_start()
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
    pub fn alternate_bases(&self) -> AlternateBases<'_> {
        self.0.alternate_bases()
    }

    /// Returns the quality score.
    pub fn quality_score(&self) -> Option<io::Result<f32>> {
        self.0.quality_score()
    }

    /// Returns the filters.
    pub fn filters(&self) -> Filters<'_> {
        self.0.filters()
    }

    /// Returns the info.
    pub fn info(&self) -> Info<'_> {
        self.0.info()
    }

    /// Returns the samples.
    pub fn samples(&self) -> Samples<'_> {
        self.0.samples()
    }
}

impl fmt::Debug for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Record")
            .field("reference_sequence_name", &self.reference_sequence_name())
            .field("variant_start", &self.variant_start())
            .field("ids", &self.ids())
            .field("reference_bases", &self.reference_bases())
            .field("alternate_bases", &self.alternate_bases())
            .field("quality_score", &self.quality_score())
            .field("filters", &self.filters())
            .field("info", &self.info())
            .field("samples", &self.samples())
            .finish()
    }
}

impl TryFrom<&[u8]> for Record {
    type Error = io::Error;

    fn try_from(buf: &[u8]) -> Result<Self, Self::Error> {
        use crate::io::Reader;

        let mut reader = Reader::new(buf);
        let mut record = Self::default();
        reader.read_record(&mut record)?;
        Ok(record)
    }
}

impl crate::variant::Record for Record {
    fn reference_sequence_name<'a, 'h: 'a>(&'a self, _: &'h Header) -> io::Result<&'a str> {
        Ok(self.reference_sequence_name())
    }

    fn variant_start(&self) -> Option<io::Result<Position>> {
        self.variant_start()
    }

    fn ids(&self) -> Box<dyn crate::variant::record::Ids + '_> {
        Box::new(self.ids())
    }

    fn reference_bases(&self) -> Box<dyn crate::variant::record::ReferenceBases + '_> {
        Box::new(self.reference_bases())
    }

    fn alternate_bases(&self) -> Box<dyn crate::variant::record::AlternateBases + '_> {
        Box::new(self.alternate_bases())
    }

    fn quality_score(&self) -> Option<io::Result<f32>> {
        self.quality_score()
    }

    fn filters(&self) -> Box<dyn crate::variant::record::Filters + '_> {
        Box::new(self.filters())
    }

    fn info(&self) -> Box<dyn crate::variant::record::Info + '_> {
        Box::new(self.info())
    }

    fn samples(&self) -> io::Result<Box<dyn crate::variant::record::Samples + '_>> {
        Ok(Box::new(self.samples()))
    }
}
