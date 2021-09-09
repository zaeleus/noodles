use std::{
    fmt,
    ops::{Bound, RangeBounds},
};

/// An interval.
pub type Interval = (Bound<i32>, Bound<i32>);

/// A mapped region.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Mapped {
    name: String,
    start: Bound<i32>,
    end: Bound<i32>,
}

impl Mapped {
    /// Creates a mapped region.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Region;
    /// assert!(Region::mapped("sq0", 5..=8).as_mapped().is_some());
    /// ```
    pub fn new<S, B>(name: S, interval: B) -> Self
    where
        S: Into<String>,
        B: RangeBounds<i32>,
    {
        Self {
            name: name.into(),
            start: interval.start_bound().cloned(),
            end: interval.end_bound().cloned(),
        }
    }

    /// Returns the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Region;
    ///
    /// assert_eq!(
    ///     Region::mapped("sq0", 5..=8).as_mapped().map(|r| r.name()),
    ///     Some("sq0")
    /// );
    /// ```
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the start position of the region (1-based).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::ops::Bound;
    /// use noodles_core::Region;
    ///
    /// assert_eq!(
    ///     Region::mapped("sq0", 5..=8).as_mapped().map(|r| r.start()),
    ///     Some(Bound::Included(5))
    /// );
    /// ```
    pub fn start(&self) -> Bound<i32> {
        self.start
    }

    /// Returns the end position of the region (1-based).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::ops::Bound;
    /// use noodles_core::Region;
    ///
    /// assert_eq!(
    ///     Region::mapped("sq0", 5..=8).as_mapped().map(|r| r.end()),
    ///     Some(Bound::Included(8))
    /// );
    /// ```
    pub fn end(&self) -> Bound<i32> {
        self.end
    }

    /// Returns the start and end positions as an interval.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::ops::Bound;
    /// use noodles_core::Region;
    ///
    /// assert_eq!(
    ///     Region::mapped("sq0", 5..=8).as_mapped().map(|r| r.interval()),
    ///     Some((Bound::Included(5), Bound::Included(8)))
    /// );
    /// ```
    pub fn interval(&self) -> Interval {
        (self.start, self.end)
    }
}

impl fmt::Display for Mapped {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.interval() {
            (Bound::Unbounded, Bound::Unbounded) => write!(f, "{}", self.name()),
            (Bound::Included(s), Bound::Unbounded) => write!(f, "{}:{}", self.name(), s),
            (Bound::Included(s), Bound::Included(e)) => write!(f, "{}:{}-{}", self.name(), s, e),
            _ => todo!(),
        }
    }
}
