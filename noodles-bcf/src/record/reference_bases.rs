/// BCF record reference bases.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceBases<'a>(&'a [u8]);

impl<'a> ReferenceBases<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }
}

impl<'a> AsRef<[u8]> for ReferenceBases<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}
