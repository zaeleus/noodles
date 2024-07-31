//! Feature record other fields.

mod value;

use bstr::BStr;

pub use self::value::Value;

/// Feature record other fields.
pub trait OtherFields {
    /// Returns whether there are any other fields.
    fn is_empty(&self) -> bool;

    /// Return the number of other fields.
    fn len(&self) -> usize;

    /// Returns an iterator over other fields.
    fn iter(&self) -> Box<dyn Iterator<Item = &BStr> + '_>;
}
