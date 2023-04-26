//! Inner VCF header other map value.

mod tag;

pub use self::tag::Tag;

use std::fmt;

use self::tag::StandardTag;
use super::{builder, Fields, Inner, Map, OtherFields, TryFromFieldsError};

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
        let mut other_fields = OtherFields::new();

        for (key, value) in fields {
            match Tag::from(key) {
                tag::ID => return Err(TryFromFieldsError::DuplicateTag),
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

        let expected = Map::<Other>::builder()
            .insert("noodles".parse()?, "vcf")
            .build()?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
