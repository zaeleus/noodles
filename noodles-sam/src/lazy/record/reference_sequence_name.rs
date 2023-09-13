/// A raw SAM record reference sequence name.
#[derive(Debug, Eq, PartialEq)]
pub struct ReferenceSequenceName<'a>(&'a [u8]);

impl<'a> ReferenceSequenceName<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }
}

impl<'a> AsRef<[u8]> for ReferenceSequenceName<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}
