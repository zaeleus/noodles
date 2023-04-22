use std::fmt;

/// A FASTQ record.
#[derive(Clone, Default, Debug, Eq, PartialEq)]
pub struct Record {
    name: Vec<u8>,
    sequence: Vec<u8>,
    description: Vec<u8>,
    quality_scores: Vec<u8>,
}

impl Record {
    /// Creates a FASTQ record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::Record;
    /// let record = Record::new("r0", "AGCT", "NDLS");
    /// ```
    pub fn new<S, T, U>(name: S, sequence: T, quality_scores: U) -> Self
    where
        S: Into<Vec<u8>>,
        T: Into<Vec<u8>>,
        U: Into<Vec<u8>>,
    {
        Self {
            name: name.into(),
            sequence: sequence.into(),
            description: Vec::new(),
            quality_scores: quality_scores.into(),
        }
    }

    /// Returns the name of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::Record;
    /// let record = Record::new("r0", "AGCT", "NDLS");
    /// assert_eq!(record.name(), b"r0");
    /// ```
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    /// Returns a mutable reference to the name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::Record;
    /// let mut record = Record::new("r0", "AGCT", "NDLS");
    /// *record.name_mut() = b"r1".to_vec();
    /// assert_eq!(record.name(), b"r1");
    /// ```
    pub fn name_mut(&mut self) -> &mut Vec<u8> {
        &mut self.name
    }

    /// Returns the sequence of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::Record;
    /// let record = Record::new("r0", "AGCT", "NDLS");
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
    /// use noodles_fastq::Record;
    /// let mut record = Record::new("r0", "AGCT", "NDLS");
    /// record.sequence_mut()[0] = b'C';
    /// assert_eq!(record.sequence(), b"CGCT");
    /// ```
    pub fn sequence_mut(&mut self) -> &mut Vec<u8> {
        &mut self.sequence
    }

    /// Returns the description of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::Record;
    /// let record = Record::new("r0", "AGCT", "NDLS");
    /// assert!(record.description().is_empty());
    /// ```
    pub fn description(&self) -> &[u8] {
        &self.description
    }

    /// Returns a mutable reference to the description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::Record;
    /// let mut record = Record::new("r0", "AGCT", "NDLS");
    /// *record.description_mut() = b"LN=4".to_vec();
    /// assert_eq!(record.description(), b"LN=4");
    /// ```
    pub fn description_mut(&mut self) -> &mut Vec<u8> {
        &mut self.description
    }

    /// Returns the quality scores of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::Record;
    /// let record = Record::new("r0", "AGCT", "NDLS");
    /// assert_eq!(record.quality_scores(), b"NDLS");
    /// ```
    pub fn quality_scores(&self) -> &[u8] {
        &self.quality_scores
    }

    /// Returns a mutable reference to the quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::Record;
    /// let mut record = Record::new("r0", "AGCT", "NDLS");
    /// *record.quality_scores_mut() = b"!!!!".to_vec();
    /// assert_eq!(record.quality_scores(), b"!!!!");
    /// ```
    pub fn quality_scores_mut(&mut self) -> &mut Vec<u8> {
        &mut self.quality_scores
    }

    // Truncates all field buffers to 0.
    pub(crate) fn clear(&mut self) {
        self.name.clear();
        self.sequence.clear();
        self.description.clear();
        self.quality_scores.clear();
    }
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("@")?;

        for &b in self.name() {
            write!(f, "{}", b as char)?;
        }

        if !self.description().is_empty() {
            write!(f, " ")?;

            for &b in self.description() {
                write!(f, "{}", b as char)?;
            }
        }

        writeln!(f)?;

        for &b in self.sequence() {
            write!(f, "{}", b as char)?;
        }

        writeln!(f)?;

        writeln!(f, "+")?;

        for &b in self.quality_scores() {
            write!(f, "{}", b as char)?;
        }

        writeln!(f)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let mut record = Record::new("r0", "ATCG", "NDLS");
        assert_eq!(record.to_string(), "@r0\nATCG\n+\nNDLS\n");

        record.description_mut().extend_from_slice(b"LN:4");
        assert_eq!(record.to_string(), "@r0 LN:4\nATCG\n+\nNDLS\n");
    }

    #[test]
    fn test_clear() {
        let mut record = Record::new("r0", "AGCT", "NDLS");
        record.clear();

        assert!(record.name().is_empty());
        assert!(record.sequence().is_empty());
        assert!(record.description().is_empty());
        assert!(record.quality_scores().is_empty());
    }
}
