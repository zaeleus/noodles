//! Inner VCF header other map value.

pub(crate) mod tag;

pub use self::tag::Tag;

use std::fmt;

use self::tag::StandardTag;
use super::{builder, Inner, Map};

/// An inner VCF header other map value.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Other;

impl Inner for Other {
    type StandardTag = StandardTag;
    type Builder = builder::Identity;
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
