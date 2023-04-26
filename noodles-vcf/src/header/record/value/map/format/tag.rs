use crate::header::record::value::map;

pub(super) type StandardTag = map::tag::TypedDescribedIndexed;

/// A VCF header format map tag.
pub type Tag = map::tag::Tag<StandardTag>;
