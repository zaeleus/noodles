//! Inner VCF header alternative allele map value.

mod builder;
pub(crate) mod tag;

pub use self::tag::Tag;

use super::{Described, Inner, Map, OtherFields};

/// An inner VCF header alternative allele map value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlternativeAllele {
    pub(crate) description: String,
}

impl Inner for AlternativeAllele {
    type StandardTag = tag::Standard;
    type Builder = builder::Builder;
}

impl Described for AlternativeAllele {
    fn description(&self) -> &str {
        &self.description
    }

    fn description_mut(&mut self) -> &mut String {
        &mut self.description
    }
}

impl Map<AlternativeAllele> {
    /// Creates a VCF header alternative allele map value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::AlternativeAllele, Map};
    /// let map = Map::<AlternativeAllele>::new("Deletion");
    /// ```
    pub fn new<D>(description: D) -> Self
    where
        D: Into<String>,
    {
        Self {
            inner: AlternativeAllele {
                description: description.into(),
            },
            other_fields: OtherFields::new(),
        }
    }
}
