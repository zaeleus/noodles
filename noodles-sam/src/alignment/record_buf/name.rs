/// An alignment record name buffer.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Name(Vec<u8>);

impl AsRef<[u8]> for Name {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<Vec<u8>> for Name {
    fn as_mut(&mut self) -> &mut Vec<u8> {
        &mut self.0
    }
}

impl From<&[u8]> for Name {
    fn from(buf: &[u8]) -> Self {
        Self::from(Vec::from(buf))
    }
}

impl<const N: usize> From<&[u8; N]> for Name {
    fn from(buf: &[u8; N]) -> Self {
        Self::from(buf.as_slice())
    }
}

impl From<Vec<u8>> for Name {
    fn from(buf: Vec<u8>) -> Self {
        Self(buf)
    }
}

impl From<Name> for Vec<u8> {
    fn from(name: Name) -> Self {
        name.0
    }
}

impl crate::alignment::record::field::Name for &Name {
    fn as_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}
