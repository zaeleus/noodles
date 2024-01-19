/// A raw SAM record name.
#[derive(Debug, Eq, PartialEq)]
pub struct Name<'a>(&'a [u8]);

impl<'a> Name<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }
}

impl<'a> crate::alignment::record::field::Name for Name<'a> {
    fn as_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl<'a> AsRef<[u8]> for Name<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}
