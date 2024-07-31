use bstr::{BStr, BString};

/// An alignment record name buffer.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Name(BString);

impl AsRef<BStr> for Name {
    fn as_ref(&self) -> &BStr {
        self.0.as_ref()
    }
}

impl AsMut<BString> for Name {
    fn as_mut(&mut self) -> &mut BString {
        &mut self.0
    }
}

impl From<&[u8]> for Name {
    fn from(buf: &[u8]) -> Self {
        Self::from(BString::from(buf))
    }
}

impl<const N: usize> From<&[u8; N]> for Name {
    fn from(buf: &[u8; N]) -> Self {
        Self::from(BString::from(buf.as_slice()))
    }
}

impl From<BString> for Name {
    fn from(buf: BString) -> Self {
        Self(buf)
    }
}

impl From<Name> for BString {
    fn from(name: Name) -> Self {
        name.0
    }
}

impl crate::alignment::record::Name for &Name {
    fn as_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}
