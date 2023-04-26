use crate::header::record::value::map;

pub(super) type StandardTag = map::tag::Identity;

/// A VCF header other map tag.
pub type Tag = map::tag::Tag<StandardTag>;
