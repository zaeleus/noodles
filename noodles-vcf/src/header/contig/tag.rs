use std::{borrow::Borrow, fmt};

/// VCF header contig record ID tag.
pub const ID: Tag = Tag::Standard(Standard::Id);

/// VCF header contig record length tag.
pub const LENGTH: Tag = Tag::Standard(Standard::Length);

/// VCF header contig record IDX tag.
pub const IDX: Tag = Tag::Standard(Standard::Idx);

/// A standard VCF header contig tag.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Standard {
    Id,
    Length,
    Idx,
}

impl Standard {
    fn new(s: &str) -> Option<Self> {
        match s {
            "ID" => Some(Self::Id),
            "length" => Some(Self::Length),
            "IDX" => Some(Self::Idx),
            _ => None,
        }
    }
}

impl fmt::Display for Standard {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Id => "ID".fmt(f),
            Self::Length => "length".fmt(f),
            Self::Idx => "IDX".fmt(f),
        }
    }
}

/// A nonstandard VCF header contig tag.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Other(String);

impl Borrow<str> for Other {
    fn borrow(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for Other {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// A VCF header contig record tag.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum Tag {
    /// A standard tag.
    Standard(Standard),
    /// A nonstandard tag.
    Other(Other),
}

impl Tag {
    /// Creates a nonstandard tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::contig::Tag;
    /// assert!(Tag::other("md5").is_some());
    /// assert!(Tag::other("ID").is_none());
    /// ```
    pub fn other(s: &str) -> Option<Other> {
        match Self::from(s) {
            Self::Standard(_) => None,
            Self::Other(tag) => Some(tag),
        }
    }
}

impl From<&str> for Tag {
    fn from(s: &str) -> Self {
        match Standard::new(s) {
            Some(tag) => Self::Standard(tag),
            None => Self::Other(Other(s.into())),
        }
    }
}

impl From<String> for Tag {
    fn from(s: String) -> Self {
        match Standard::new(&s) {
            Some(tag) => Self::Standard(tag),
            None => Self::Other(Other(s)),
        }
    }
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Standard(inner) => inner.fmt(f),
            Self::Other(inner) => inner.fmt(f),
        }
    }
}
