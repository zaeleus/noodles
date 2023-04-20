//! Inner VCF header other map value.

use std::fmt;

use super::{builder, tag, Fields, Inner, Map, TryFromFieldsError};

type StandardTag = tag::Identity;

/// A VCF header other map tag.
pub type Tag = tag::Tag<StandardTag>;

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

impl TryFrom<Fields> for Map<Other> {
    type Error = TryFromFieldsError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut other_fields = super::init_other_fields();

        for (key, value) in fields {
            match Tag::from(key) {
                Tag::Standard(StandardTag::Id) => return Err(TryFromFieldsError::DuplicateTag),
                Tag::Other(t) => super::insert_other_field(&mut other_fields, t, value)?,
            }
        }

        Ok(Self {
            inner: Other,
            other_fields,
        })
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

    #[test]
    fn test_try_from_fields_for_map_other() -> Result<(), Box<dyn std::error::Error>> {
        let actual = Map::<Other>::try_from(vec![(String::from("noodles"), String::from("vcf"))])?;

        let noodles_tag = match Tag::from(String::from("noodles")) {
            Tag::Other(tag) => tag,
            Tag::Standard(_) => panic!("invalid tag"),
        };

        let expected = Map::<Other>::builder().insert(noodles_tag, "vcf").build()?;
        assert_eq!(actual, expected);

        Ok(())
    }
}
