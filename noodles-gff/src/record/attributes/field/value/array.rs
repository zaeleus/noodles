use std::{borrow::Cow, fmt, io};

use crate::record::attributes::field::percent_decode;

/// A raw GFF record attributes field array value.
#[derive(Eq, PartialEq)]
pub struct Array<'a>(&'a str);

impl<'a> Array<'a> {
    pub(super) fn new(s: &'a str) -> Self {
        Self(s)
    }

    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<Cow<'a, str>>> {
        const DELIMITER: char = ',';

        self.0
            .split(DELIMITER)
            .map(|s| percent_decode(s).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
    }
}

impl AsRef<str> for Array<'_> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

impl fmt::Debug for Array<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl<'a> crate::feature::record::attributes::field::value::Array<'a> for Array<'a> {
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Cow<'a, str>>> + 'a> {
        Box::new(self.iter())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iter() -> io::Result<()> {
        let array = Array::new("8,13%2C21");
        let actual: Vec<_> = array.iter().collect::<io::Result<_>>()?;
        assert_eq!(actual, [Cow::from("8"), Cow::from("13,21")]);
        Ok(())
    }
}
