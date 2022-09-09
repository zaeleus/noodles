use std::fmt;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub struct ContentId(i32);

impl fmt::Display for ContentId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl From<i32> for ContentId {
    fn from(n: i32) -> Self {
        Self(n)
    }
}

impl From<ContentId> for i32 {
    fn from(content_id: ContentId) -> Self {
        content_id.0
    }
}
