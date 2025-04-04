//! GFF record attributes.

pub mod field;

use std::{borrow::Cow, io};

use bstr::BStr;
use indexmap::IndexMap;

use self::field::{Tag, Value};

/// GFF record attributes.
///
/// Attributes are extra data attached to a GFF record. They are represented as a typed map, where
/// each key ([`Tag`]) is associated with a typed [`Value`].
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Attributes(IndexMap<Tag, Value>);

impl Attributes {
    /// Returns whether there are any entries.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of entries.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns the value at the given tag.
    pub fn get(&self, tag: &[u8]) -> Option<&Value> {
        self.0.get(tag)
    }
}

impl AsRef<IndexMap<Tag, Value>> for Attributes {
    fn as_ref(&self) -> &IndexMap<Tag, Value> {
        &self.0
    }
}

impl AsMut<IndexMap<Tag, Value>> for Attributes {
    fn as_mut(&mut self) -> &mut IndexMap<Tag, Value> {
        &mut self.0
    }
}

impl Extend<(Tag, Value)> for Attributes {
    fn extend<T: IntoIterator<Item = (Tag, Value)>>(&mut self, iter: T) {
        self.0.extend(iter);
    }
}

impl FromIterator<(Tag, Value)> for Attributes {
    fn from_iter<T: IntoIterator<Item = (Tag, Value)>>(iter: T) -> Self {
        let mut attributes = Self::default();
        attributes.extend(iter);
        attributes
    }
}

impl crate::feature::record::Attributes for Attributes {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn get(
        &self,
        tag: &[u8],
    ) -> Option<io::Result<crate::feature::record::attributes::field::Value<'_>>> {
        self.get(tag).map(|value| Ok(value.into()))
    }

    fn iter(
        &self,
    ) -> Box<
        dyn Iterator<
                Item = io::Result<(
                    Cow<'_, BStr>,
                    crate::feature::record::attributes::field::Value<'_>,
                )>,
            > + '_,
    > {
        Box::new(
            self.0
                .iter()
                .map(|(tag, value)| Ok((Cow::from(tag), value.into()))),
        )
    }
}

impl crate::feature::record::Attributes for &Attributes {
    fn is_empty(&self) -> bool {
        Attributes::is_empty(self)
    }

    fn get(
        &self,
        tag: &[u8],
    ) -> Option<io::Result<crate::feature::record::attributes::field::Value<'_>>> {
        Attributes::get(self, tag).map(|value| Ok(value.into()))
    }

    fn iter(
        &self,
    ) -> Box<
        dyn Iterator<
                Item = io::Result<(
                    Cow<'_, BStr>,
                    crate::feature::record::attributes::field::Value<'_>,
                )>,
            > + '_,
    > {
        Box::new(
            self.0
                .iter()
                .map(|(tag, value)| Ok((Cow::from(tag), value.into()))),
        )
    }
}
