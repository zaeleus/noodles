/// BCF record genotypes.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Genotypes {
    buf: Vec<u8>,
    format_count: usize,
    sample_count: usize,
}

impl Genotypes {
    /// Returns the number of samples.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Genotypes;
    /// let genotypes = Genotypes::default();
    /// assert_eq!(genotypes.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.sample_count
    }

    /// Returns whether there are any samples.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Genotypes;
    /// let genotypes = Genotypes::default();
    /// assert!(genotypes.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the number of fields per sample.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Genotypes;
    /// let genotypes = Genotypes::default();
    /// assert_eq!(genotypes.format_count(), 0);
    /// ```
    pub fn format_count(&self) -> usize {
        self.format_count
    }

    pub(crate) fn set_format_count(&mut self, format_count: usize) {
        self.format_count = format_count;
    }

    pub(crate) fn set_sample_count(&mut self, sample_count: usize) {
        self.sample_count = sample_count;
    }
}

impl AsRef<[u8]> for Genotypes {
    fn as_ref(&self) -> &[u8] {
        &self.buf
    }
}

impl AsMut<Vec<u8>> for Genotypes {
    fn as_mut(&mut self) -> &mut Vec<u8> {
        &mut self.buf
    }
}
