/// A gzip index.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Index(Vec<(u64, u64)>);

impl AsRef<[(u64, u64)]> for Index {
    fn as_ref(&self) -> &[(u64, u64)] {
        &self.0
    }
}

impl From<Vec<(u64, u64)>> for Index {
    fn from(index: Vec<(u64, u64)>) -> Self {
        Self(index)
    }
}
