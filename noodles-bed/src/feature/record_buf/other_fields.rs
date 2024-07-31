//! Feature record other fields buffer.

mod value;

pub use self::value::Value;

/// A feature record other fields buffer.
#[derive(Clone, Default, Debug, PartialEq)]
pub struct OtherFields(Vec<Value>);

impl AsRef<[Value]> for OtherFields {
    fn as_ref(&self) -> &[Value] {
        &self.0
    }
}

impl AsMut<Vec<Value>> for OtherFields {
    fn as_mut(&mut self) -> &mut Vec<Value> {
        &mut self.0
    }
}

impl crate::feature::record::OtherFields for OtherFields {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(
        &self,
    ) -> Box<dyn Iterator<Item = crate::feature::record::other_fields::Value<'_>> + '_> {
        Box::new(self.0.iter().map(|value| value.into()))
    }
}

impl crate::feature::record::OtherFields for &OtherFields {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(
        &self,
    ) -> Box<dyn Iterator<Item = crate::feature::record::other_fields::Value<'_>> + '_> {
        Box::new(self.0.iter().map(|value| value.into()))
    }
}
