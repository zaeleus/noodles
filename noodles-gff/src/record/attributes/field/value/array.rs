use std::fmt;

/// A raw GFF record attributes field array value.
#[derive(Eq, PartialEq)]
pub struct Array<'a>(&'a str);

impl<'a> Array<'a> {
    pub(super) fn new(s: &'a str) -> Self {
        Self(s)
    }

    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = &'a str> {
        const DELIMITER: char = ',';
        self.0.split(DELIMITER)
    }
}

impl<'a> AsRef<str> for Array<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

impl<'a> fmt::Debug for Array<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iter() {
        let array = Array::new("nd,ls");
        let actual: Vec<_> = array.iter().collect();
        assert_eq!(actual, ["nd", "ls"]);
    }
}
