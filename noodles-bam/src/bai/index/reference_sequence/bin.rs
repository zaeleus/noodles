//! BAM index bin and fields.

mod chunk;

pub use self::chunk::Chunk;

/// A bin in a BAM index reference sequence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Bin {
    bin: u32,
    chunks: Vec<Chunk>,
}

impl Bin {
    /// Creates a new bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, Vec::new());
    /// ```
    pub fn new(bin: u32, chunks: Vec<Chunk>) -> Self {
        Self { bin, chunks }
    }

    /// Returns the bin index of this bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, Vec::new());
    /// assert_eq!(bin.bin(), 10946);
    /// ```
    pub fn bin(&self) -> u32 {
        self.bin
    }

    /// Returns the list of chunks in the bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, Vec::new());
    /// assert!(bin.chunks().is_empty());
    /// ```
    pub fn chunks(&self) -> &[Chunk] {
        &self.chunks
    }
}
