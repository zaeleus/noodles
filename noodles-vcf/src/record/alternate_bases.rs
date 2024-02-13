use std::fmt;

/// VCF record alternate bases.
pub struct AlternateBases<'a>(&'a str);

const DELIMITER: char = ',';

impl<'a> AlternateBases<'a> {
    pub(super) fn new(src: &'a str) -> Self {
        Self(src)
    }

    /// Returns whether there are any alternate bases.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of alternate bases.
    pub fn len(&self) -> usize {
        if self.is_empty() {
            0
        } else {
            self.iter().count()
        }
    }

    /// Returns an iterator over alternate bases.
    pub fn iter(&self) -> impl Iterator<Item = &'a str> + '_ {
        self.0.split(DELIMITER)
    }
}

impl<'a> fmt::Debug for AlternateBases<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl<'a> AsRef<str> for AlternateBases<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}
