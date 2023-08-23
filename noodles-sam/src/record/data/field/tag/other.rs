use super::LENGTH;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[doc(hidden)]
pub struct Other(pub(super) [u8; LENGTH]);

impl AsRef<[u8; LENGTH]> for Other {
    fn as_ref(&self) -> &[u8; LENGTH] {
        &self.0
    }
}
