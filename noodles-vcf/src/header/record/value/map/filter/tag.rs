use crate::header::record::value::map;

pub(crate) type Standard = map::tag::DescribedIndexed;

/// A VCF header filter map tag.
pub type Tag = map::tag::Tag<Standard>;

// For some reason, using the `Tag` type alias produces a `nontrivial_structural_match` warning
// when pattern matching, so it's avoided here.
pub(crate) const ID: Tag = map::tag::Tag::Standard(Standard::Id);
pub(crate) const DESCRIPTION: Tag = map::tag::Tag::Standard(Standard::Description);
pub(crate) const IDX: Tag = map::tag::Tag::Standard(Standard::Idx);
