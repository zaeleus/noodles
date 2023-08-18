//! Inner VCF header filter map value.

pub(crate) mod tag;

pub use self::tag::Tag;

use std::fmt;

use indexmap::IndexMap;

use super::{builder, Described, Indexed, Inner, Map};

/// An inner VCF header filter map value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Filter {
    pub(crate) description: String,
    pub(crate) idx: Option<usize>,
}

impl Inner for Filter {
    type StandardTag = tag::Standard;
    type Builder = builder::DescribedIndexed;
}

impl Described for Filter {
    fn description(&self) -> &str {
        &self.description
    }

    fn description_mut(&mut self) -> &mut String {
        &mut self.description
    }
}

impl Indexed for Filter {
    fn idx(&self) -> Option<usize> {
        self.idx
    }

    fn idx_mut(&mut self) -> &mut Option<usize> {
        &mut self.idx
    }
}

impl Map<Filter> {
    /// Creates a default VCF header filter map value for "PASS".
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Filter, Map};
    /// let actual = Map::<Filter>::pass();
    /// let expected = Map::<Filter>::new("All filters passed");
    /// assert_eq!(actual, expected);
    /// ```
    pub fn pass() -> Self {
        Self {
            inner: Filter {
                description: String::from("All filters passed"),
                idx: None,
            },
            other_fields: IndexMap::new(),
        }
    }

    /// Creates a VCF header filter map value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Filter, Map};
    /// let map = Map::<Filter>::new("Quality below 10");
    /// ```
    pub fn new<D>(description: D) -> Self
    where
        D: Into<String>,
    {
        Self {
            inner: Filter {
                description: description.into(),
                idx: None,
            },
            other_fields: IndexMap::new(),
        }
    }
}

impl fmt::Display for Map<Filter> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        super::fmt_display_description_field(f, self.description())?;
        super::fmt_display_other_fields(f, self.other_fields())?;

        if let Some(idx) = self.idx() {
            super::fmt_display_idx_field(f, idx)?;
        }

        Ok(())
    }
}

impl builder::Inner<Filter> for builder::DescribedIndexed {
    fn build(self) -> Result<Filter, builder::BuildError> {
        let description = self
            .description
            .ok_or(builder::BuildError::MissingField("Description"))?;

        Ok(Filter {
            description,
            idx: self.idx,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let map = Map::<Filter>::pass();
        let expected = r#",Description="All filters passed""#;
        assert_eq!(map.to_string(), expected);
    }
}
