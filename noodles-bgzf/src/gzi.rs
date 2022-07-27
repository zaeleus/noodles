//! A module for the [GZI Index]. A GZI index contains pairs of compressed and uncompressed offsets
//! in a BGZF file. Values in the index are stored as little-endian 64-bit unsigned integers.
//!
//! [GZI Index]: http://www.htslib.org/doc/bgzip.html#GZI_FORMAT

pub use self::reader::Reader;

mod reader;

/// A GZI Index contains a number of entries, representing pairs of compressed and uncompressed
/// offsets in a BGZF file.
#[derive(Debug, PartialEq)]
pub struct Index {
    number_entries: u64,
    offsets: Vec<(u64, u64)>,
}

impl Index {
    /// Creates a GZI Index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::gzi;
    ///
    /// let index = gzi::Index::new(vec![(4668, 21294)]);
    /// ```
    pub fn new(offsets: Vec<(u64, u64)>) -> Self {
        Self {
            number_entries: offsets.len() as u64,
            offsets,
        }
    }

    /// Returns the number of entries.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::gzi;
    ///
    /// let index = gzi::Index::new(vec![(4668, 21294)]);
    /// assert_eq!(index.number_entries(), 1);
    /// ```
    pub fn number_entries(&self) -> u64 {
        self.number_entries
    }

    /// Returns the compressed and uncompressed offset pairs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::gzi;
    ///
    /// let index = gzi::Index::new(vec![(4668, 21294)]);
    /// assert_eq!(index.offsets(), [(4668, 21294)]);
    /// ```
    pub fn offsets(&self) -> &[(u64, u64)] {
        &self.offsets
    }
}

#[cfg(test)]
mod tests {
    use crate::gzi::Index;

    #[test]
    fn test_number_entries() {
        let index = Index::new(vec![(4668, 21294)]);
        assert_eq!(index.number_entries(), 1);
    }

    #[test]
    fn test_offsets() {
        let index = Index::new(vec![(4668, 21294)]);
        assert_eq!(index.offsets(), [(4668, 21294)]);
    }
}
