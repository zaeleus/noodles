//! Inner VCF header other map value.

pub(crate) mod tag;

pub use self::tag::Tag;

use std::fmt;

use super::{builder, Inner, Map};

/// An inner VCF header other map value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Other {
    pub(crate) id_tag: Tag,
}

impl Inner for Other {
    type StandardTag = tag::Standard;
    type Builder = builder::Identity;
}

impl Default for Other {
    fn default() -> Self {
        Self { id_tag: tag::ID }
    }
}

impl Map<Other> {
    /// Creates a nonstandard VCF header map value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Other, Map};
    /// let map = Map::<Other>::new();
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    pub(crate) fn id_tag(&self) -> &Tag {
        &self.inner.id_tag
    }
}

impl fmt::Display for Map<Other> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        super::fmt_display_other_fields(f, self.other_fields())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let map = Map::<Other>::new();
        assert!(map.to_string().is_empty());
    }
}
