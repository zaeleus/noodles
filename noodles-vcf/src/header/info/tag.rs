use std::{borrow::Borrow, fmt};

/* /// A VCF header INFO record ID tag.
pub const ID: Tag = Tag::Standard(Standard::Id);

/// A VCF header INFO record number tag.
pub const NUMBER: Tag = Tag::Standard(Standard::Number);

/// A VCF header INFO record type tag.
pub const TYPE: Tag = Tag::Standard(Standard::Type);

/// A VCF header INFO record description tag.
pub const DESCRIPTION: Tag = Tag::Standard(Standard::Description);

/// A VCF header INFO record IDX tag.
pub const IDX: Tag = Tag::Standard(Standard::Idx); */

/// A standard VCF header INFO tag.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Standard {
    Id,
    Number,
    Type,
    Description,
    Idx,
}

impl Standard {
    fn new(s: &str) -> Option<Self> {
        match s {
            "ID" => Some(Self::Id),
            "Number" => Some(Self::Number),
            "Type" => Some(Self::Type),
            "Description" => Some(Self::Description),
            "IDX" => Some(Self::Idx),
            _ => None,
        }
    }
}

impl AsRef<str> for Standard {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => "ID",
            Self::Number => "Number",
            Self::Type => "Type",
            Self::Description => "Description",
            Self::Idx => "IDX",
        }
    }
}

impl fmt::Display for Standard {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.as_ref().fmt(f)
    }
}

/// A nonstandard VCF header INFO tag.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Other(String);

impl AsRef<str> for Other {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl Borrow<str> for Other {
    fn borrow(&self) -> &str {
        self.as_ref()
    }
}

impl fmt::Display for Other {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// A VCF header INFO record tag.
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

impl AsRef<str> for Tag {
    fn as_ref(&self) -> &str {
        match self {
            Self::Standard(tag) => tag.as_ref(),
            Self::Other(tag) => tag.as_ref(),
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
