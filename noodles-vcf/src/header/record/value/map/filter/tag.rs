use crate::header::record::value::map;

pub(super) type StandardTag = map::tag::DescribedIndexed;

/// A VCF header filter map tag.
pub type Tag = map::tag::Tag<StandardTag>;
