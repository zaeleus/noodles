use bytes::Bytes;

/// A FASTA record sequence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Sequence(Bytes);

impl Sequence {
    /// Returns the length of the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::record::Sequence;
    /// let sequence = Sequence::default();
    /// assert_eq!(sequence.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns whether the sequence is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::record::Sequence;
    /// let sequence = Sequence::default();
    /// assert!(sequence.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl AsRef<[u8]> for Sequence {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl From<Vec<u8>> for Sequence {
    fn from(data: Vec<u8>) -> Self {
        Self(Bytes::from(data))
    }
}
