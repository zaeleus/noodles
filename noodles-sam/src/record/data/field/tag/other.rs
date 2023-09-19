use super::LENGTH;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[doc(hidden)]
pub struct Other(pub(super) [u8; LENGTH]);

impl Other {
    pub(super) const fn new(buf: [u8; LENGTH]) -> Option<Self> {
        if buf[0].is_ascii_alphabetic() && buf[1].is_ascii_alphanumeric() {
            Some(Self(buf))
        } else {
            None
        }
    }
}

impl AsRef<[u8; LENGTH]> for Other {
    fn as_ref(&self) -> &[u8; LENGTH] {
        &self.0
    }
}
