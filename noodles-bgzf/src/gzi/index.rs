use std::io;

use crate::VirtualPosition;

/// A gzip index (GZI).
///
/// A gzip index holds compressed-uncompressed position pairs.
///
/// Like the physical index, this does _not_ include the position of the first block, which is
/// implicitly at 0.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Index(Vec<(u64, u64)>);

impl Index {
    /// Returns the virtual position at the given uncompressed position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::{self as bgzf, gzi};
    ///
    /// let index = gzi::Index::default();
    /// assert_eq!(index.query(0)?, bgzf::VirtualPosition::default());
    ///
    /// let index = gzi::Index::from(vec![(8, 21), (13, 55)]);
    /// assert_eq!(index.query(0)?, bgzf::VirtualPosition::default());
    /// assert_eq!(index.query(13)?, bgzf::VirtualPosition::try_from((0, 13))?);
    /// assert_eq!(index.query(34)?, bgzf::VirtualPosition::try_from((8, 13))?);
    /// assert_eq!(index.query(89)?, bgzf::VirtualPosition::try_from((13, 34))?);
    /// Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn query(&self, pos: u64) -> io::Result<VirtualPosition> {
        let i = self.0.partition_point(|r| r.1 <= pos);

        let (compressed_pos, uncompressed_pos) = if i == 0 { (0, 0) } else { self.0[i - 1] };

        let block_data_pos = u16::try_from(pos - uncompressed_pos)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        VirtualPosition::try_from((compressed_pos, block_data_pos))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
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
