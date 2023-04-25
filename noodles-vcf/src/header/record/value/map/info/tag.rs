use crate::header::record::value::map;

pub(super) type StandardTag = map::tag::TypedDescribedIndexed;

/// A VCF header info map tag.
pub type Tag = map::tag::Tag<StandardTag>;
