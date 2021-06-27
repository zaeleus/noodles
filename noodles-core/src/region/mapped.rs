use std::{
    fmt,
    ops::{Bound, RangeBounds},
};

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
            start: bound_cloned(interval.start_bound()),
            end: bound_cloned(interval.end_bound()),
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
    pub fn interval(&self) -> (Bound<i32>, Bound<i32>) {
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

// TODO: https://github.com/rust-lang/rust/issues/61356
fn bound_cloned<T>(bound: Bound<&T>) -> Bound<T>
where
    T: Clone,
{
    match bound {
        Bound::Included(v) => Bound::Included(v.clone()),
        Bound::Excluded(v) => Bound::Excluded(v.clone()),
        Bound::Unbounded => Bound::Unbounded,
    }
}
