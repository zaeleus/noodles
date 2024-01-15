//! FASTA record and definition.

pub mod definition;
pub mod sequence;

pub use self::{definition::Definition, sequence::Sequence};

/// A FASTA record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Record {
    definition: Definition,
    sequence: Sequence,
}

impl Record {
    /// Creates a FASTA record from a definition and sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::{self as fasta, record::{Definition, Sequence}};
    /// let definition = Definition::new("sq0", None);
    /// let sequence = Sequence::from(b"ACGT".to_vec());
    /// let record = fasta::Record::new(definition, sequence);
    /// ```
    pub fn new(definition: Definition, sequence: Sequence) -> Self {
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
    /// use noodles_fasta::{self as fasta, record::{Definition, Sequence}};
    ///
    /// let definition = Definition::new("sq0", None);
    /// let sequence = Sequence::from(b"ACGT".to_vec());
    /// let record = fasta::Record::new(definition.clone(), sequence);
    ///
    /// assert_eq!(record.definition(), &definition);
    /// ```
    pub fn definition(&self) -> &Definition {
        &self.definition
    }

    /// Returns the record name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::{self as fasta, record::{Definition, Sequence}};
    ///
    /// let definition = Definition::new("sq0", None);
    /// let sequence = Sequence::from(b"ACGT".to_vec());
    /// let record = fasta::Record::new(definition, sequence);
    ///
    /// assert_eq!(record.name(), b"sq0");
    /// ```
    pub fn name(&self) -> &[u8] {
        self.definition.name()
    }

    /// Returns the description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::{self as fasta, record::{Definition, Sequence}};
    ///
    /// let definition = Definition::new("sq0", Some(Vec::from("LN:4")));
    /// let sequence = Sequence::from(b"ACGT".to_vec());
    /// let record = fasta::Record::new(definition, sequence);
    ///
    /// assert_eq!(record.description(), Some(&b"LN:4"[..]));
    /// ```
    pub fn description(&self) -> Option<&[u8]> {
        self.definition.description()
    }

    /// Returns the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::{self as fasta, record::{Definition, Sequence}};
    ///
    /// let definition = Definition::new("sq0", None);
    /// let sequence = Sequence::from(b"ACGT".to_vec());
    /// let record = fasta::Record::new(definition, sequence.clone());
    ///
    /// assert_eq!(record.sequence(), &sequence);
    /// ```
    pub fn sequence(&self) -> &Sequence {
        &self.sequence
    }
}
