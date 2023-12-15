/// A raw SAM record read name.
#[derive(Debug, Eq, PartialEq)]
pub struct ReadName<'a>(&'a [u8]);

impl<'a> ReadName<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }
}

impl<'a> AsRef<[u8]> for ReadName<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}
