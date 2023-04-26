use crate::header::record::value::map;

pub(super) type StandardTag = map::tag::DescribedIndexed;

/// A VCF header filter map tag.
pub type Tag = map::tag::Tag<StandardTag>;

// For some reason, using the `Tag` type alias produces a `nontrivial_structural_match` warning
// when pattern matching, so it's avoided here.
pub(crate) const ID: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Id);
pub(super) const DESCRIPTION: Tag =
    map::tag::Tag::<StandardTag>::Standard(StandardTag::Description);
pub(super) const IDX: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Idx);
