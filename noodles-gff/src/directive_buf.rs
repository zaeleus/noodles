//! GFF directives.

pub mod key;
pub mod value;

use crate::Directive;

pub use self::value::Value;

/// A GFF directive.
///
/// This is also called a pragma or metadata.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DirectiveBuf {
    key: String,
    value: Option<Value>,
}

impl DirectiveBuf {
    /// Creates a directive buffer.
    pub fn new<K>(key: K, value: Option<Value>) -> Self
    where
        K: Into<String>,
    {
        Self {
            key: key.into(),
            value,
        }
    }

    /// Returns the key.
    pub fn key(&self) -> &str {
        &self.key
    }

    /// Returns the value.
    pub fn value(&self) -> Option<&Value> {
        self.value.as_ref()
    }
}

impl From<Directive<'_>> for DirectiveBuf {
    fn from(directive: Directive<'_>) -> Self {
        Self::new(directive.key(), directive.value().map(Value::from))
    }
}
