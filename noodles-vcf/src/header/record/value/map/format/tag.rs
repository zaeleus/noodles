use crate::header::record::value::map;

pub(crate) type StandardTag = map::tag::TypedDescribedIndexed;

/// A VCF header format map tag.
pub type Tag = map::tag::Tag<StandardTag>;

// For some reason, using the `Tag` type alias produces a `nontrivial_structural_match` warning
// when pattern matching, so it's avoided here.
pub(crate) const ID: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Id);
pub(crate) const NUMBER: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Number);
pub(crate) const TYPE: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Type);
pub(crate) const DESCRIPTION: Tag =
    map::tag::Tag::<StandardTag>::Standard(StandardTag::Description);
pub(crate) const IDX: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Idx);
