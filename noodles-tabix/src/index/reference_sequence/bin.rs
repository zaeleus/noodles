//! Tabix index bin and fields.

mod builder;
mod chunk;

pub use self::chunk::Chunk;

pub(crate) use self::builder::Builder;

// MAX_BIN (2019-04-09)
pub(crate) const MAX_ID: usize = ((1 << 18) - 1) / 7 + 1;

/// A tabix index reference sequence bin.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Bin {
    id: u32,
    chunks: Vec<Chunk>,
}

impl Bin {
    pub(crate) fn builder() -> Builder {
        Builder::default()
    }

    /// Creates a new bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, Vec::new());
    /// ```
    pub fn new(id: u32, chunks: Vec<Chunk>) -> Self {
        Self { id, chunks }
    }

    /// Returns the bin ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, Vec::new());
    /// assert_eq!(bin.id(), 10946);
    /// ```
    pub fn id(&self) -> u32 {
        self.id
    }

    /// Returns the list of chunks in the bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, Vec::new());
    /// assert!(bin.chunks().is_empty());
    /// ```
    pub fn chunks(&self) -> &[Chunk] {
        &self.chunks
    }
}
