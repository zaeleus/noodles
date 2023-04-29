use super::LENGTH;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[doc(hidden)]
pub struct Other(pub(super) [u8; LENGTH]);
