//! FASTA record and definition.

pub mod definition;

pub use self::definition::Definition;

/// A FASTA record.
#[derive(Clone, Debug, Eq, PartialEq)]
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

    /// Returns the record definition.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    ///
    /// let definition = fasta::record::Definition::new(String::from("sq0"), None);
    /// let sequence = b"ACGT".to_vec();
    /// let record = fasta::Record::new(definition.clone(), sequence);
    ///
    /// assert_eq!(record.definition(), &definition);
    /// ```
    pub fn definition(&self) -> &Definition {
        &self.definition
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
    /// assert_eq!(record.reference_sequence_name(), "sq0");
    /// ```
    #[deprecated(since = "0.2.0", note = "Use `name` instead.")]
    pub fn reference_sequence_name(&self) -> &str {
        self.definition.name()
    }

    /// Returns the record name.
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
        self.definition.name()
    }

    /// Returns the description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    ///
    /// let definition = fasta::record::Definition::new(
    ///     String::from("sq0"),
    ///     Some(String::from("LN:4"))
    /// );
    ///
    /// let sequence = b"ACGT".to_vec();
    /// let record = fasta::Record::new(definition, sequence);
    ///
    /// assert_eq!(record.description(), Some("LN:4"));
    /// ```
    pub fn description(&self) -> Option<&str> {
        self.definition.description()
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
