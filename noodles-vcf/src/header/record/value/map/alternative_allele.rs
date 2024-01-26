//! Inner VCF header alternative allele map value.

mod builder;
pub(crate) mod tag;

pub use self::tag::Tag;

use std::{error, fmt};

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
    /// use noodles_vcf::{
    ///     header::record::value::{map::AlternativeAllele, Map},
    ///     record::alternate_bases::allele::{
    ///         symbol::{structural_variant::Type, StructuralVariant},
    ///         Symbol,
    ///     },
    /// };
    ///
    /// let id = Symbol::StructuralVariant(StructuralVariant::from(Type::Deletion));
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

/// An error returned when a raw ALT record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is missing.
    MissingField(Tag),
    /// A tag is duplicated.
    DuplicateTag(Tag),
    /// The ID is invalid.
    InvalidId(crate::record::alternate_bases::allele::symbol::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidId(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(tag) => write!(f, "missing field: {tag}"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
            Self::InvalidId(_) => write!(f, "invalid ID"),
        }
    }
}
