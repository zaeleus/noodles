//! GFF record attributes field value.

use std::{
    error, fmt, iter, mem,
    str::{self, FromStr},
};

const DELIMITER: char = ',';

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
    /// use noodles_gff::record::attributes::field::Value;
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
    /// use noodles_gff::record::attributes::field::Value;
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
    /// use noodles_gff::record::attributes::field::Value;
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

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use super::percent_encode;

        for (i, value) in self.iter().enumerate() {
            if i > 0 {
                DELIMITER.fmt(f)?;
            }

            percent_encode(value).fmt(f)?;
        }

        Ok(())
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

/// An error returned when a raw GFF record attribute field value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid(str::Utf8Error),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(_) => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Value {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some((a, b)) = s.split_once(DELIMITER) {
            iter::once(a)
                .chain(b.split(DELIMITER))
                .map(decode_value)
                .collect::<Result<_, _>>()
                .map(Self::Array)
        } else {
            decode_value(s).map(Self::String)
        }
    }
}

fn decode_value(s: &str) -> Result<String, ParseError> {
    use super::percent_decode;

    percent_decode(s)
        .map(|t| t.into_owned())
        .map_err(ParseError::Invalid)
}

impl<'a> IntoIterator for &'a Value {
    type Item = &'a String;
    type IntoIter = Box<dyn Iterator<Item = Self::Item> + 'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let value = Value::from("gene0");
        assert_eq!(value.to_string(), "gene0");

        let value = Value::from("13,21");
        assert_eq!(value.to_string(), "13%2C21");

        let value = Value::from(vec![String::from("13"), String::from("21")]);
        assert_eq!(value.to_string(), "13,21");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("gene0".parse(), Ok(Value::from("gene0")));
        assert_eq!("13%2C21".parse(), Ok(Value::from("13,21")));
        assert_eq!(
            "13,21".parse(),
            Ok(Value::from(vec![String::from("13"), String::from("21")]))
        );
    }
}
