//! FASTQ record.

mod definition;

pub use self::definition::Definition;

use std::fmt;

/// A FASTQ record.
#[derive(Clone, Default, Eq, PartialEq)]
pub struct Record {
    definition: Definition,
    sequence: Vec<u8>,
    quality_scores: Vec<u8>,
}

impl Record {
    /// Creates a FASTQ record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::{self as fastq, record::Definition};
    /// let record = fastq::Record::new(Definition::new("r0", ""), "AGCT", "NDLS");
    /// ```
    pub fn new<S, Q>(definition: Definition, sequence: S, quality_scores: Q) -> Self
    where
        S: Into<Vec<u8>>,
        Q: Into<Vec<u8>>,
    {
        Self {
            definition,
            sequence: sequence.into(),
            quality_scores: quality_scores.into(),
        }
    }

    /// Returns the record definition.
    pub fn definition(&self) -> &Definition {
        &self.definition
    }

    /// Returns a mutable reference to the record definition.
    pub fn definition_mut(&mut self) -> &mut Definition {
        &mut self.definition
    }

    /// Returns the name of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::{self as fastq, record::Definition};
    /// let record = fastq::Record::new(Definition::new("r0", ""), "AGCT", "NDLS");
    /// assert_eq!(record.name(), b"r0");
    /// ```
    pub fn name(&self) -> &[u8] {
        self.definition.name()
    }

    /// Returns a mutable reference to the name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::{self as fastq, record::Definition};
    /// let mut record = fastq::Record::new(Definition::new("r0", ""), "AGCT", "NDLS");
    /// *record.name_mut() = b"r1".to_vec();
    /// assert_eq!(record.name(), b"r1");
    /// ```
    pub fn name_mut(&mut self) -> &mut Vec<u8> {
        self.definition.name_mut()
    }

    /// Returns the description of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::{self as fastq, record::Definition};
    /// let record = fastq::Record::new(Definition::new("r0", ""), "AGCT", "NDLS");
    /// assert!(record.description().is_empty());
    /// ```
    pub fn description(&self) -> &[u8] {
        self.definition.description()
    }

    /// Returns a mutable reference to the description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::{self as fastq, record::Definition};
    /// let mut record = fastq::Record::new(Definition::new("r0", ""), "AGCT", "NDLS");
    /// *record.description_mut() = b"LN=4".to_vec();
    /// assert_eq!(record.description(), b"LN=4");
    /// ```
    pub fn description_mut(&mut self) -> &mut Vec<u8> {
        self.definition.description_mut()
    }

    /// Returns the sequence of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::{self as fastq, record::Definition};
    /// let record = fastq::Record::new(Definition::new("r0", ""), "AGCT", "NDLS");
    /// assert_eq!(record.sequence(), b"AGCT");
    /// ```
    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    /// Returns a mutable reference to the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::{self as fastq, record::Definition};
    /// let mut record = fastq::Record::new(Definition::new("r0", ""), "AGCT", "NDLS");
    /// record.sequence_mut()[0] = b'C';
    /// assert_eq!(record.sequence(), b"CGCT");
    /// ```
    pub fn sequence_mut(&mut self) -> &mut Vec<u8> {
        &mut self.sequence
    }

    /// Returns the quality scores of the record.
    ///
    /// The encoding of these scores are considered to be unknown; and it is up to the caller to
    /// decode them, if necessary.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::{self as fastq, record::Definition};
    /// let record = fastq::Record::new(Definition::new("r0", ""), "AGCT", "NDLS");
    /// assert_eq!(record.quality_scores(), b"NDLS");
    /// ```
    pub fn quality_scores(&self) -> &[u8] {
        &self.quality_scores
    }

    /// Returns a mutable reference to the quality scores.
    ///
    /// The encoding of these scores are considered to be unknown; and it is up to the caller to
    /// encode/decode them, if necessary.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::{self as fastq, record::Definition};
    /// let mut record = fastq::Record::new(Definition::new("r0", ""), "AGCT", "NDLS");
    /// *record.quality_scores_mut() = b"!!!!".to_vec();
    /// assert_eq!(record.quality_scores(), b"!!!!");
    /// ```
    pub fn quality_scores_mut(&mut self) -> &mut Vec<u8> {
        &mut self.quality_scores
    }

    // Truncates all field buffers to 0.
    pub(crate) fn clear(&mut self) {
        self.definition.clear();
        self.sequence.clear();
        self.quality_scores.clear();
    }
}

impl fmt::Debug for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use std::str;

        let name = str::from_utf8(self.name());
        let description = str::from_utf8(self.description());
        let sequence = str::from_utf8(self.sequence());
        let quality_scores = str::from_utf8(self.quality_scores());

        f.debug_struct("Record")
            .field("name", &name)
            .field("description", &description)
            .field("sequence", &sequence)
            .field("quality_scores", &quality_scores)
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clear() {
        let mut record = Record::new(Definition::new("r0", ""), "AGCT", "NDLS");
        record.clear();

        assert!(record.name().is_empty());
        assert!(record.description().is_empty());
        assert!(record.sequence().is_empty());
        assert!(record.quality_scores().is_empty());
    }
}
