use crate::header::record::value::map;

pub(super) type StandardTag = map::tag::Described;

/// A VCF header alternative allele map tag.
pub type Tag = map::tag::Tag<StandardTag>;
