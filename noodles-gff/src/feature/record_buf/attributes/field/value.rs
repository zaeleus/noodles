//! GFF record attributes field value.

use std::{borrow::Cow, io, iter};

use bstr::{BStr, BString};

/// A GFF record attribute field value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Value {
    /// A string.
    String(BString),
    /// An array.
    Array(Vec<BString>),
}

impl Value {
    /// Returns the value as a string, if the value is a string.
    ///
    /// # Examples
    ///
    /// ```
    /// use bstr::{BStr, BString};
    /// use noodles_gff::feature::record_buf::attributes::field::Value;
    ///
    /// let value = Value::from("ndls");
    /// assert_eq!(value.as_string(), Some(BStr::new("ndls")));
    ///
    /// let value = Value::from(vec![BString::from("ndls0"), BString::from("ndls1")]);
    /// assert!(value.as_string().is_none());
    /// ```
    pub fn as_string(&self) -> Option<&BStr> {
        match self {
            Self::String(value) => Some(value.as_ref()),
            Self::Array(_) => None,
        }
    }

    /// Returns the value as an array, if the value is an array.
    ///
    /// # Examples
    ///
    /// ```
    /// use bstr::BString;
    /// use noodles_gff::feature::record_buf::attributes::field::Value;
    ///
    /// let raw_values = vec![BString::from("ndls0"), BString::from("ndls1")];
    /// let value = Value::from(raw_values.clone());
    /// assert_eq!(value.as_array(), Some(&raw_values[..]));
    ///
    /// let value = Value::from("ndls");
    /// assert!(value.as_array().is_none());
    /// ```
    pub fn as_array(&self) -> Option<&[BString]> {
        match self {
            Self::String(_) => None,
            Self::Array(values) => Some(values),
        }
    }

    /// An iterator over values.
    ///
    /// # Examples
    ///
    /// ```
    /// use bstr::BString;
    /// use noodles_gff::feature::record_buf::attributes::field::Value;
    ///
    /// let value = Value::from("ndls");
    /// let mut iter = value.iter();
    /// assert_eq!(iter.next(), Some(&BString::from("ndls")));
    /// assert!(iter.next().is_none());
    ///
    /// let value = Value::from(vec![BString::from("ndls0"), BString::from("ndls1")]);
    /// let mut iter = value.iter();
    /// assert_eq!(iter.next(), Some(&BString::from("ndls0")));
    /// assert_eq!(iter.next(), Some(&BString::from("ndls1")));
    /// assert!(iter.next().is_none());
    /// ```
    pub fn iter(&self) -> Box<dyn Iterator<Item = &BString> + '_> {
        match self {
            Self::String(value) => Box::new(iter::once(value)),
            Self::Array(values) => Box::new(values.iter()),
        }
    }
}

impl Extend<BString> for Value {
    fn extend<T: IntoIterator<Item = BString>>(&mut self, iter: T) {
        match self {
            Self::String(value) => {
                let mut values = vec![value.clone()];
                values.extend(iter);
                *self = Self::Array(values);
            }
            Self::Array(values) => values.extend(iter),
        }
    }
}

impl From<&str> for Value {
    fn from(s: &str) -> Self {
        Self::String(s.into())
    }
}

impl From<String> for Value {
    fn from(s: String) -> Self {
        Self::String(s.into())
    }
}

impl From<Vec<BString>> for Value {
    fn from(values: Vec<BString>) -> Self {
        Self::Array(values)
    }
}

impl<'a> IntoIterator for &'a Value {
    type Item = &'a BString;
    type IntoIter = Box<dyn Iterator<Item = Self::Item> + 'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a> From<&'a Value> for crate::feature::record::attributes::field::Value<'a> {
    fn from(value_buf: &'a Value) -> Self {
        match value_buf {
            Value::String(value) => Self::String(Cow::from(value)),
            Value::Array(values) => Self::Array(Box::new(Array(values))),
        }
    }
}

struct Array<'a>(&'a [BString]);

impl<'a> crate::feature::record::attributes::field::value::Array<'a> for Array<'a> {
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Cow<'a, BStr>>> + 'a> {
        Box::new(self.0.iter().map(|value| Ok(Cow::from(value))))
    }
}
