/// A FASTQ record.
#[derive(Default, Debug, Eq, PartialEq)]
pub struct Record {
    name: Vec<u8>,
    sequence: Vec<u8>,
    plus_line: Vec<u8>,
    quality_scores: Vec<u8>,
}

impl Record {
    /// Creates a new FASTQ record with the given lines.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::Record;
    ///
    /// let record = Record::new("@noodles/1", "AGCT", "+", "NDLS");
    ///
    /// assert_eq!(record.name(), b"@noodles/1");
    /// assert_eq!(record.sequence(), b"AGCT");
    /// assert_eq!(record.plus_line(), b"+");
    /// assert_eq!(record.quality_scores(), b"NDLS");
    /// ```
    pub fn new<S, T, U, V>(name: S, sequence: T, plus_line: U, quality_scores: V) -> Self
    where
        S: Into<Vec<u8>>,
        T: Into<Vec<u8>>,
        U: Into<Vec<u8>>,
        V: Into<Vec<u8>>,
    {
        Self {
            name: name.into(),
            sequence: sequence.into(),
            plus_line: plus_line.into(),
            quality_scores: quality_scores.into(),
        }
    }

    pub fn name(&self) -> &[u8] {
        &self.name
    }

    pub fn name_mut(&mut self) -> &mut Vec<u8> {
        &mut self.name
    }

    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    pub fn sequence_mut(&mut self) -> &mut Vec<u8> {
        &mut self.sequence
    }

    pub fn plus_line(&self) -> &[u8] {
        &self.plus_line
    }

    pub fn plus_line_mut(&mut self) -> &mut Vec<u8> {
        &mut self.plus_line
    }

    pub fn quality_scores(&self) -> &[u8] {
        &self.quality_scores
    }

    pub fn quality_scores_mut(&mut self) -> &mut Vec<u8> {
        &mut self.quality_scores
    }

    /// Truncates all line buffers to 0.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::Record;
    ///
    /// let mut record = Record::new("@noodles/1", "AGCT", "+", "NDLS");
    /// record.clear();
    ///
    /// assert!(record.name().is_empty());
    /// assert!(record.sequence().is_empty());
    /// assert!(record.plus_line().is_empty());
    /// assert!(record.quality_scores().is_empty());
    /// ```
    pub fn clear(&mut self) {
        self.name.clear();
        self.sequence.clear();
        self.plus_line.clear();
        self.quality_scores.clear();
    }
}
