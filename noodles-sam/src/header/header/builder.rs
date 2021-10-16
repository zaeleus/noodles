use std::collections::HashMap;

use super::{GroupOrder, Header, SortOrder, SubsortOrder, Tag, Version};

/// A SAM header header builder.
#[derive(Debug, Default)]
pub struct Builder {
    version: Option<Version>,
    sort_order: Option<SortOrder>,
    group_order: Option<GroupOrder>,
    subsort_order: Option<SubsortOrder>,
    fields: HashMap<Tag, String>,
}

impl Builder {
    /// Sets a format version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::{Header, Version};
    /// let header = Header::builder().set_version(Version::new(1, 6)).build();
    /// assert_eq!(header.version(), Version::new(1, 6));
    /// ```
    pub fn set_version(mut self, version: Version) -> Self {
        self.version = Some(version);
        self
    }

    /// Sets a sort order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::{Header, SortOrder};
    /// let header = Header::builder().set_sort_order(SortOrder::Coordinate).build();
    /// assert_eq!(header.sort_order(), Some(SortOrder::Coordinate));
    /// ```
    pub fn set_sort_order(mut self, sort_order: SortOrder) -> Self {
        self.sort_order = Some(sort_order);
        self
    }

    /// Sets a group order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::{GroupOrder, Header};
    /// let header = Header::builder().set_group_order(GroupOrder::Reference).build();
    /// assert_eq!(header.group_order(), Some(GroupOrder::Reference));
    /// ```
    pub fn set_group_order(mut self, group_order: GroupOrder) -> Self {
        self.group_order = Some(group_order);
        self
    }

    /// Sets a subsort order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::{Header, SubsortOrder};
    ///
    /// let header = Header::builder()
    ///     .set_subsort_order(SubsortOrder::Coordinate(vec![String::from("MI")]))
    ///     .build();
    ///
    /// assert_eq!(header.subsort_order(), Some(&SubsortOrder::Coordinate(vec![String::from("MI")])));
    /// ```
    pub fn set_subsort_order(mut self, subsort_order: SubsortOrder) -> Self {
        self.subsort_order = Some(subsort_order);
        self
    }

    /// Inserts a tag-raw value pair.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::{Header, Tag};
    ///
    /// let zn = Tag::Other([b'z', b'n']);
    ///
    /// let header = Header::builder()
    ///     .insert(zn.clone(), String::from("noodles"))
    ///     .build();
    ///
    /// assert_eq!(header.fields().get(&zn), Some(&String::from("noodles")));
    /// ```
    pub fn insert<I>(mut self, tag: Tag, value: I) -> Self
    where
        I: Into<String>,
    {
        self.fields.insert(tag, value.into());
        self
    }

    /// Builds a header header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::Header;
    /// let header = Header::builder().build();
    /// ```
    pub fn build(self) -> Header {
        Header {
            version: self.version.unwrap_or_default(),
            sort_order: self.sort_order,
            group_order: self.group_order,
            subsort_order: self.subsort_order,
            fields: self.fields,
        }
    }
}
