/// A gzip index.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Index(Vec<(u64, u64)>);

impl Index {
    /// Returns the record of the block that contains the given position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::gzi;
    ///
    /// let index = gzi::Index::default();
    /// assert!(index.query(0).is_none());
    ///
    /// let index = gzi::Index::from(vec![(0, 0), (8, 21), (13, 55)]);
    /// assert_eq!(index.query(0), Some((0, 0)));
    /// assert_eq!(index.query(13), Some((0, 0)));
    /// assert_eq!(index.query(34), Some((8, 21)));
    /// assert_eq!(index.query(89), Some((13, 55)));
    /// ```
    pub fn query(&self, pos: u64) -> Option<(u64, u64)> {
        if self.0.is_empty() {
            None
        } else {
            let i = self.0.partition_point(|r| r.1 <= pos);
            // SAFETY: `i` is > 0.
            Some(self.0[i - 1])
        }
    }
}

impl AsRef<[(u64, u64)]> for Index {
    fn as_ref(&self) -> &[(u64, u64)] {
        &self.0
    }
}

impl From<Vec<(u64, u64)>> for Index {
    fn from(index: Vec<(u64, u64)>) -> Self {
        Self(index)
    }
}
