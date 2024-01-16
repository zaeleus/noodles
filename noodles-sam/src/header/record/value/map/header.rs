//! SAM header record header map value.

mod builder;
pub mod group_order;
pub mod sort_order;
pub mod subsort_order;
pub mod tag;
pub mod version;

pub use self::{
    group_order::GroupOrder, sort_order::SortOrder, subsort_order::SubsortOrder, tag::Tag,
    version::Version,
};

use self::builder::Builder;
use super::{Inner, Map, OtherFields};

/// A SAM header record header map value.
///
/// The header describes file-level metadata. The format version is guaranteed to be set.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Header {
    pub(crate) version: Version,
}

impl Inner for Header {
    type StandardTag = tag::Standard;
    type Builder = Builder;
}

impl Map<Header> {
    /// Creates a SAM header record header map value with a format version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{
    ///     map::{self, header::Version},
    ///     Map,
    /// };
    ///
    /// let header = Map::<map::Header>::new(Version::new(1, 6));
    /// ```
    pub fn new(version: Version) -> Self {
        Self {
            inner: Header { version },
            other_fields: OtherFields::new(),
        }
    }

    /// Returns the format version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{
    ///     map::{self, header::Version},
    ///     Map,
    /// };
    ///
    /// let header = Map::<map::Header>::new(Version::new(1, 6));
    /// assert_eq!(header.version(), Version::new(1, 6));
    /// ```
    pub fn version(&self) -> Version {
        self.inner.version
    }

    /// Returns a mutable reference to the format version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::{self, header::Version}, Map};
    /// let mut header = Map::<map::Header>::default();
    /// *header.version_mut() = Version::new(1, 5);
    /// assert_eq!(header.version(), Version::new(1, 5));
    /// ```
    pub fn version_mut(&mut self) -> &mut Version {
        &mut self.inner.version
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let header = Map::<Header>::default();
        assert_eq!(header.version(), Version::default());
    }
}
