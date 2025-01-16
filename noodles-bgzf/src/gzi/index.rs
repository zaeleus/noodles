/// A gzip index (GZI).
///
/// A gzip index holds compressed-uncompressed position pairs.
///
/// Like this physical index, this does _not_ include the position of the first block, which is
/// implicity at 0.
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
    /// assert_eq!(index.query(0), (0, 0));
    ///
    /// let index = gzi::Index::from(vec![(8, 21), (13, 55)]);
    /// assert_eq!(index.query(0), (0, 0));
    /// assert_eq!(index.query(13), (0, 0));
    /// assert_eq!(index.query(34), (8, 21));
    /// assert_eq!(index.query(89), (13, 55));
    /// ```
    pub fn query(&self, pos: u64) -> (u64, u64) {
        let i = self.0.partition_point(|r| r.1 <= pos);

        if i == 0 {
            (0, 0)
        } else {
            self.0[i - 1]
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
