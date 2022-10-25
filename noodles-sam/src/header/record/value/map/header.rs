//! SAM header record header map value.

mod builder;
pub mod group_order;
pub mod sort_order;
pub mod subsort_order;
mod tag;
pub mod version;

pub use self::{
    group_order::GroupOrder, sort_order::SortOrder, subsort_order::SubsortOrder, version::Version,
};

use std::fmt;

use self::builder::Builder;
use super::{Fields, Inner, Map, OtherFields, TryFromFieldsError};

type StandardTag = tag::Standard;
type Tag = super::tag::Tag<StandardTag>;

/// A SAM header record header map value.
///
/// The header describes file-level metadata. The format version is guaranteed to be set.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Header {
    version: Version,
    sort_order: Option<SortOrder>,
    group_order: Option<GroupOrder>,
    subsort_order: Option<SubsortOrder>,
}

impl Inner for Header {
    type StandardTag = StandardTag;
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
            inner: Header {
                version,
                ..Default::default()
            },
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

    /// Returns the sort order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map, Map};
    /// let header = Map::<map::Header>::default();
    /// assert!(header.sort_order().is_none());
    /// ```
    pub fn sort_order(&self) -> Option<SortOrder> {
        self.inner.sort_order
    }

    /// Returns a mutable reference to the sort order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::{self, header::SortOrder}, Map};
    /// let mut header = Map::<map::Header>::default();
    /// *header.sort_order_mut() = Some(SortOrder::Coordinate);
    /// assert_eq!(header.sort_order(), Some(SortOrder::Coordinate));
    /// ```
    pub fn sort_order_mut(&mut self) -> &mut Option<SortOrder> {
        &mut self.inner.sort_order
    }

    /// Returns the group order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map, Map};
    /// let header = Map::<map::Header>::default();
    /// assert!(header.group_order().is_none());
    /// ```
    pub fn group_order(&self) -> Option<GroupOrder> {
        self.inner.group_order
    }

    /// Returns a mutable reference to the group order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::{self, header::GroupOrder}, Map};
    /// let mut header = Map::<map::Header>::default();
    /// *header.group_order_mut() = Some(GroupOrder::None);
    /// assert_eq!(header.group_order(), Some(GroupOrder::None));
    /// ```
    pub fn group_order_mut(&mut self) -> &mut Option<GroupOrder> {
        &mut self.inner.group_order
    }

    /// Returns the subsort order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map, Map};
    /// let header = Map::<map::Header>::default();
    /// assert!(header.subsort_order().is_none());
    /// ```
    pub fn subsort_order(&self) -> Option<&SubsortOrder> {
        self.inner.subsort_order.as_ref()
    }

    /// Returns a mutable reference to the subsort order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::{self, header::SubsortOrder}, Map};
    /// let subsort_order: SubsortOrder = "coordinate:queryname".parse()?;
    /// let mut header = Map::<map::Header>::default();
    /// *header.subsort_order_mut() = Some(subsort_order.clone());
    /// assert_eq!(header.subsort_order(), Some(&subsort_order));
    /// # Ok::<_, noodles_sam::header::record::value::map::header::subsort_order::ParseError>(())
    /// ```
    pub fn subsort_order_mut(&mut self) -> &mut Option<SubsortOrder> {
        &mut self.inner.subsort_order
    }
}

impl fmt::Display for Map<Header> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "VN:{}", self.version())?;

        if let Some(sort_order) = self.sort_order() {
            write!(f, "\tSO:{}", sort_order)?;
        }

        if let Some(group_order) = self.group_order() {
            write!(f, "\tGO:{}", group_order)?;
        }

        if let Some(subsort_order) = self.subsort_order() {
            write!(f, "\tSS:{}", subsort_order)?;
        }

        super::fmt_display_other_fields(f, self.other_fields())?;

        Ok(())
    }
}

impl TryFrom<Fields> for Map<Header> {
    type Error = TryFromFieldsError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut other_fields = super::init_other_fields(fields.len());

        let mut version = None;
        let mut sort_order = None;
        let mut group_order = None;
        let mut subsort_order = None;

        for (key, value) in fields {
            let tag = key.parse().map_err(|_| TryFromFieldsError::InvalidTag)?;

            match tag {
                Tag::Standard(StandardTag::Version) => {
                    version = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("VN"))?;
                }
                Tag::Standard(StandardTag::SortOrder) => {
                    sort_order = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("SO"))?;
                }
                Tag::Standard(StandardTag::GroupOrder) => {
                    group_order = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("GO"))?;
                }
                Tag::Standard(StandardTag::SubsortOrder) => {
                    subsort_order = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("SS"))?;
                }
                Tag::Other(t) => super::insert_other_field(&mut other_fields, t, value)?,
            }
        }

        let version = version.ok_or(TryFromFieldsError::MissingField("VN"))?;

        Ok(Self {
            inner: Header {
                version,
                sort_order,
                group_order,
                subsort_order,
            },
            other_fields,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::header::record::value::map::builder::BuildError;

    #[test]
    fn test_default() {
        let header = Map::<Header>::default();
        assert_eq!(header.version(), Version::default());
        assert!(header.sort_order().is_none());
        assert!(header.group_order().is_none());
        assert!(header.subsort_order().is_none());
        assert!(header.other_fields().is_empty());
    }

    #[test]
    fn test_fmt() -> Result<(), BuildError> {
        let header = Map::<Header>::builder()
            .set_version(Version::new(1, 6))
            .set_sort_order(SortOrder::Unsorted)
            .set_group_order(GroupOrder::Query)
            .build()?;

        assert_eq!(header.to_string(), "VN:1.6\tSO:unsorted\tGO:query");

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_header_with_missing_version() -> Result<(), TryFromFieldsError>
    {
        let fields = vec![(String::from("SO"), String::from("coordinate"))];

        assert_eq!(
            Map::<Header>::try_from(fields),
            Err(TryFromFieldsError::MissingField("VN"))
        );

        Ok(())
    }
}
