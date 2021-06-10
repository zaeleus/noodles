//! SAM header record value.

use std::{error, fmt};

use indexmap::IndexMap;

/// An ordered map of raw tag-value pairs.
pub type Fields = IndexMap<String, String>;

/// A SAM header record value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Value {
    /// A string.
    String(String),
    /// A map of fields.
    Map(Fields),
}

/// An error returned when raw SAM header record fields fail to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromIteratorError {
    /// The input is empty.
    Empty,
    /// A tag is duplicated.
    DuplicateTag(String),
}

impl error::Error for TryFromIteratorError {}

impl fmt::Display for TryFromIteratorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {}", tag),
        }
    }
}

impl Value {
    /// Performs a conversion from a `(String, String)` iterator to a `Value::Map`.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::{value, Value};
    ///
    /// let actual = Value::try_from_iter(vec![("SN", "sq0"), ("LN", "8")])?;
    /// let expected = Value::Map(vec![
    ///     (String::from("SN"), String::from("sq0")),
    ///     (String::from("LN"), String::from("8")),
    /// ].into_iter().collect());
    ///
    /// assert_eq!(actual, expected);
    /// # Ok::<(), value::TryFromIteratorError>(())
    /// ```
    pub fn try_from_iter<I, K, V>(iter: I) -> Result<Self, TryFromIteratorError>
    where
        I: IntoIterator<Item = (K, V)>,
        K: Into<String>,
        V: Into<String>,
    {
        let mut fields = Fields::new();

        for (k, v) in iter {
            let key = k.into();
            let value = v.into();

            if fields.insert(key.clone(), value).is_some() {
                return Err(TryFromIteratorError::DuplicateTag(key));
            }
        }

        if fields.is_empty() {
            return Err(TryFromIteratorError::Empty);
        }

        Ok(Value::Map(fields))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_iter() {
        assert_eq!(
            Value::try_from_iter(Vec::<(String, String)>::new()),
            Err(TryFromIteratorError::Empty)
        );
    }
}
