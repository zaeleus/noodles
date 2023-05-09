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

use std::{error, fmt};

use self::{
    builder::Builder,
    tag::{StandardTag, Tag},
};
use super::{Fields, Inner, Map, OtherFields};
use crate::header::parser::Context;

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
        write!(f, "{}:{}", tag::VERSION, self.version())?;

        if let Some(sort_order) = self.sort_order() {
            write!(f, "\t{}:{sort_order}", tag::SORT_ORDER)?;
        }

        if let Some(group_order) = self.group_order() {
            write!(f, "\t{}:{group_order}", tag::GROUP_ORDER)?;
        }

        if let Some(subsort_order) = self.subsort_order() {
            write!(f, "\t{}:{subsort_order}", tag::SUBSORT_ORDER)?;
        }

        super::fmt_display_other_fields(f, self.other_fields())?;

        Ok(())
    }
}

/// An error returned when a raw header header record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is missing.
    MissingField(Tag),
    /// A tag is invalid.
    InvalidTag(super::tag::ParseError),
    /// A tag is duplicated.
    DuplicateTag(Tag),
    /// The version is invalid.
    InvalidVersion(version::ParseError),
    /// The sort order is invalid.
    InvalidSortOrder(sort_order::ParseError),
    /// The group order is invalid.
    InvalidGroupOrder(group_order::ParseError),
    /// The subsort order is invalid.
    InvalidSubsortOrder(subsort_order::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidTag(e) => Some(e),
            Self::InvalidVersion(e) => Some(e),
            Self::InvalidSortOrder(e) => Some(e),
            Self::InvalidGroupOrder(e) => Some(e),
            Self::InvalidSubsortOrder(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(tag) => write!(f, "missing field: {tag}"),
            Self::InvalidTag(_) => write!(f, "invalid tag"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
            Self::InvalidVersion(_) => write!(f, "invalid version"),
            Self::InvalidSortOrder(_) => write!(f, "invalid sort order"),
            Self::InvalidGroupOrder(_) => write!(f, "invalid group order"),
            Self::InvalidSubsortOrder(_) => write!(f, "invalid subsort order"),
        }
    }
}

impl TryFrom<Fields> for Map<Header> {
    type Error = ParseError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        Self::try_from((&Context::default(), fields))
    }
}

impl TryFrom<(&Context, Fields)> for Map<Header> {
    type Error = ParseError;

    fn try_from((ctx, fields): (&Context, Fields)) -> Result<Self, Self::Error> {
        let mut version = None;
        let mut sort_order = None;
        let mut group_order = None;
        let mut subsort_order = None;

        let mut other_fields = OtherFields::new();

        for (key, value) in fields {
            let tag = key.parse().map_err(ParseError::InvalidTag)?;

            match tag {
                tag::VERSION => parse_version(&value)
                    .and_then(|v| try_replace(&mut version, ctx, tag::VERSION, v))?,
                tag::SORT_ORDER => parse_sort_order(&value)
                    .and_then(|v| try_replace(&mut sort_order, ctx, tag::SORT_ORDER, v))?,
                tag::GROUP_ORDER => parse_group_order(&value)
                    .and_then(|v| try_replace(&mut group_order, ctx, tag::GROUP_ORDER, v))?,
                tag::SUBSORT_ORDER => parse_subsort_order(&value)
                    .and_then(|v| try_replace(&mut subsort_order, ctx, tag::SUBSORT_ORDER, v))?,
                Tag::Other(t) => try_insert(&mut other_fields, ctx, t, value)?,
            }
        }

        let version = version.ok_or(ParseError::MissingField(tag::VERSION))?;

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

fn parse_version(s: &str) -> Result<Version, ParseError> {
    s.parse().map_err(ParseError::InvalidVersion)
}

fn parse_sort_order(s: &str) -> Result<SortOrder, ParseError> {
    s.parse().map_err(ParseError::InvalidSortOrder)
}

fn parse_group_order(s: &str) -> Result<GroupOrder, ParseError> {
    s.parse().map_err(ParseError::InvalidGroupOrder)
}

fn parse_subsort_order(s: &str) -> Result<SubsortOrder, ParseError> {
    s.parse().map_err(ParseError::InvalidSubsortOrder)
}

fn try_replace<T>(
    option: &mut Option<T>,
    ctx: &Context,
    tag: Tag,
    value: T,
) -> Result<(), ParseError> {
    if option.replace(value).is_some() && !ctx.allow_duplicate_tags() {
        Err(ParseError::DuplicateTag(tag))
    } else {
        Ok(())
    }
}

fn try_insert(
    other_fields: &mut OtherFields<StandardTag>,
    ctx: &Context,
    tag: super::tag::Other<StandardTag>,
    value: String,
) -> Result<(), ParseError> {
    if other_fields.insert(tag, value).is_some() && !ctx.allow_duplicate_tags() {
        Err(ParseError::DuplicateTag(Tag::Other(tag)))
    } else {
        Ok(())
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
    fn test_try_from_fields_for_map_header_with_missing_version() {
        assert_eq!(
            Map::<Header>::try_from(vec![(String::from("SO"), String::from("coordinate"))]),
            Err(ParseError::MissingField(tag::VERSION))
        );
    }
}
