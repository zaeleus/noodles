/// A FASTQ record containing a name, sequence, plus line, and quality.
#[derive(Default, Debug)]
pub struct Record {
    name: Vec<u8>,
    sequence: Vec<u8>,
    plus_line: Vec<u8>,
    quality: Vec<u8>,
}

impl Record {
    /// Creates a new FASTQ record with the given lines.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles::formats::fastq::Record;
    ///
    /// let record = Record::new("@noodles/1", "AGCT", "+", "abcd");
    ///
    /// assert_eq!(record.name(), b"@noodles/1");
    /// assert_eq!(record.sequence(), b"AGCT");
    /// assert_eq!(record.plus_line(), b"+");
    /// assert_eq!(record.quality(), b"abcd");
    /// ```
    pub fn new<S, T, U, V>(
        name: S,
        sequence: T,
        plus_line: U,
        quality: V,
    ) -> Record
    where
        S: Into<Vec<u8>>,
        T: Into<Vec<u8>>,
        U: Into<Vec<u8>>,
        V: Into<Vec<u8>>,
    {
        Record {
            name: name.into(),
            sequence: sequence.into(),
            plus_line: plus_line.into(),
            quality: quality.into(),
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

    pub fn quality(&self) -> &[u8] {
        &self.quality
    }

    pub fn quality_mut(&mut self) -> &mut Vec<u8> {
        &mut self.quality
    }

    /// Truncates all line buffers to 0.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles::formats::fastq::Record;
    ///
    /// let mut record = Record::new("@noodles/1", "AGCT", "+", "abcd");
    /// record.clear();
    ///
    /// assert!(record.name().is_empty());
    /// assert!(record.sequence().is_empty());
    /// assert!(record.plus_line().is_empty());
    /// assert!(record.quality().is_empty());
    /// ```
    pub fn clear(&mut self) {
        self.name.clear();
        self.sequence.clear();
        self.plus_line.clear();
        self.quality.clear();
    }
}
