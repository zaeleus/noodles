//! GFF record attributes field value.

use std::{iter, mem};

/// A GFF record attribute field value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Value {
    /// A string.
    String(String),
    /// An array.
    Array(Vec<String>),
}

impl Value {
    /// Returns the value as a string, if the value is a string.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::record_buf::attributes::field::Value;
    ///
    /// let value = Value::from("ndls");
    /// assert_eq!(value.as_string(), Some("ndls"));
    ///
    /// let value = Value::from(vec![String::from("ndls0"), String::from("ndls1")]);
    /// assert!(value.as_string().is_none());
    /// ```
    pub fn as_string(&self) -> Option<&str> {
        match self {
            Self::String(value) => Some(value),
            Self::Array(_) => None,
        }
    }

    /// Returns the value as an array, if the value is an array.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::record_buf::attributes::field::Value;
    ///
    /// let raw_values = vec![String::from("ndls0"), String::from("ndls1")];
    /// let value = Value::from(raw_values.clone());
    /// assert_eq!(value.as_array(), Some(&raw_values[..]));
    ///
    /// let value = Value::from("ndls");
    /// assert!(value.as_array().is_none());
    /// ```
    pub fn as_array(&self) -> Option<&[String]> {
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
    /// use noodles_gff::record_buf::attributes::field::Value;
    ///
    /// let value = Value::from("ndls");
    /// let mut iter = value.iter();
    /// assert_eq!(iter.next(), Some(&String::from("ndls")));
    /// assert!(iter.next().is_none());
    ///
    /// let value = Value::from(vec![String::from("ndls0"), String::from("ndls1")]);
    /// let mut iter = value.iter();
    /// assert_eq!(iter.next(), Some(&String::from("ndls0")));
    /// assert_eq!(iter.next(), Some(&String::from("ndls1")));
    /// assert!(iter.next().is_none());
    /// ```
    pub fn iter(&self) -> Box<dyn Iterator<Item = &String> + '_> {
        match self {
            Self::String(value) => Box::new(iter::once(value)),
            Self::Array(values) => Box::new(values.iter()),
        }
    }
}

impl Extend<String> for Value {
    fn extend<T: IntoIterator<Item = String>>(&mut self, iter: T) {
        match self {
            Self::String(value) => {
                let mut values = vec![value.clone()];
                values.extend(iter);
                mem::swap(self, &mut Self::Array(values));
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
        Self::String(s)
    }
}

impl From<Vec<String>> for Value {
    fn from(values: Vec<String>) -> Self {
        Self::Array(values)
    }
}

impl<'a> IntoIterator for &'a Value {
    type Item = &'a String;
    type IntoIter = Box<dyn Iterator<Item = Self::Item> + 'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}
