use crate::header::record::value::map;

pub(super) type StandardTag = map::tag::TypedDescribedIndexed;

/// A VCF header info map tag.
pub type Tag = map::tag::Tag<StandardTag>;

// For some reason, using the `Tag` type alias produces a `nontrivial_structural_match` warning
// when pattern matching, so it's avoided here.
pub(crate) const ID: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Id);
pub(super) const NUMBER: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Number);
pub(super) const TYPE: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Type);
pub(super) const DESCRIPTION: Tag =
    map::tag::Tag::<StandardTag>::Standard(StandardTag::Description);
pub(super) const IDX: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Idx);
