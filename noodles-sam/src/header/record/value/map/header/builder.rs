use super::{GroupOrder, Header, SortOrder, SubsortOrder, Version};
use crate::header::record::value::map::{self, builder::BuildError};

/// A SAM header header builder.
#[derive(Debug, Default)]
pub struct Builder {
    version: Option<Version>,
    sort_order: Option<SortOrder>,
    group_order: Option<GroupOrder>,
    subsort_order: Option<SubsortOrder>,
}

impl map::Builder<Header> {
    /// Sets a format version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{
    ///     map::{self, header::Version},
    ///     Map,
    /// };
    ///
    /// let version = Version::new(1, 6);
    ///
    /// let header = Map::<map::Header>::builder()
    ///     .set_version(version.clone())
    ///     .build()?;
    ///
    /// assert_eq!(header.version(), version);
    /// Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_version(mut self, version: Version) -> Self {
        self.inner.version = Some(version);
        self
    }

    /// Sets a sort order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{
    ///     map::{self, header::SortOrder},
    ///     Map,
    /// };
    ///
    /// let header = Map::<map::Header>::builder()
    ///     .set_sort_order(SortOrder::Coordinate)
    ///     .build()?;
    ///
    /// assert_eq!(header.sort_order(), Some(SortOrder::Coordinate));
    /// Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_sort_order(mut self, sort_order: SortOrder) -> Self {
        self.inner.sort_order = Some(sort_order);
        self
    }

    /// Sets a group order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{
    ///     map::{self, header::GroupOrder},
    ///     Map,
    /// };
    ///
    /// let header = Map::<map::Header>::builder()
    ///     .set_group_order(GroupOrder::Reference)
    ///     .build()?;
    ///
    /// assert_eq!(header.group_order(), Some(GroupOrder::Reference));
    /// Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_group_order(mut self, group_order: GroupOrder) -> Self {
        self.inner.group_order = Some(group_order);
        self
    }

    /// Sets a subsort order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{
    ///     map::{self, header::SubsortOrder},
    ///     Map,
    /// };
    ///
    /// let subsort_order = SubsortOrder::Coordinate(vec![String::from("MI")]);
    ///
    /// let header = Map::<map::Header>::builder()
    ///     .set_subsort_order(subsort_order.clone())
    ///     .build()?;
    ///
    /// assert_eq!(header.subsort_order(), Some(&subsort_order));
    /// Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_subsort_order(mut self, subsort_order: SubsortOrder) -> Self {
        self.inner.subsort_order = Some(subsort_order);
        self
    }
}

impl map::builder::Inner<Header> for Builder {
    fn build(self) -> Result<Header, BuildError> {
        Ok(Header {
            version: self.version.unwrap_or_default(),
            sort_order: self.sort_order,
            group_order: self.group_order,
            subsort_order: self.subsort_order,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert!(builder.version.is_none());
        assert!(builder.sort_order.is_none());
        assert!(builder.group_order.is_none());
        assert!(builder.subsort_order.is_none());
    }
}
