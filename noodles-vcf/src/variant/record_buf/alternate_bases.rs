//! VCF record alternate bases.

/// VCF record alternate bases (`ALT`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AlternateBases(Vec<String>);

impl AlternateBases {
    /// Return whether there are any alternate alleles.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of alternate alleles.
    pub fn len(&self) -> usize {
        self.0.len()
    }
}

impl AsRef<[String]> for AlternateBases {
    fn as_ref(&self) -> &[String] {
        &self.0
    }
}

impl AsMut<Vec<String>> for AlternateBases {
    fn as_mut(&mut self) -> &mut Vec<String> {
        &mut self.0
    }
}

impl From<Vec<String>> for AlternateBases {
    fn from(alleles: Vec<String>) -> Self {
        Self(alleles)
    }
}
