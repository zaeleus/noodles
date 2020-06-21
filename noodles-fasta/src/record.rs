//! FASTA record and definition.

pub mod definition;

pub use self::definition::Definition;

/// A FASTA record.
pub struct Record {
    definition: Definition,
    sequence: Vec<u8>,
}

impl Record {
    /// Creates a FASTA record from a definition and sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// let definition = fasta::record::Definition::new(String::from("sq0"), None);
    /// let sequence = b"ACGT".to_vec();
    /// let record = fasta::Record::new(definition, sequence);
    /// ```
    pub fn new(definition: Definition, sequence: Vec<u8>) -> Self {
        Self {
            definition,
            sequence,
        }
    }

    /// Returns the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    ///
    /// let definition = fasta::record::Definition::new(String::from("sq0"), None);
    /// let sequence = b"ACGT".to_vec();
    /// let record = fasta::Record::new(definition, sequence);
    ///
    /// assert_eq!(record.name(), "sq0");
    /// ```
    pub fn name(&self) -> &str {
        &self.definition.reference_sequence_name()
    }

    /// Returns the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    ///
    /// let definition = fasta::record::Definition::new(String::from("sq0"), None);
    /// let sequence = b"ACGT".to_vec();
    /// let record = fasta::Record::new(definition, sequence);
    ///
    /// assert_eq!(record.sequence(), b"ACGT");
    /// ```
    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }
}
