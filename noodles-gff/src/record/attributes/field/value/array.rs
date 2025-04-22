use std::{borrow::Cow, fmt, io};

use bstr::{BStr, ByteSlice};

use crate::record::attributes::field::percent_decode;

/// A GFF record attributes field array value.
#[derive(Eq, PartialEq)]
pub struct Array<'a>(&'a [u8]);

impl<'a> Array<'a> {
    pub(super) fn new(s: &'a [u8]) -> Self {
        Self(s)
    }

    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = Cow<'a, BStr>> {
        const DELIMITER: u8 = b',';
        self.0.split(|b| *b == DELIMITER).map(|s| percent_decode(s))
    }
}

impl AsRef<BStr> for Array<'_> {
    fn as_ref(&self) -> &BStr {
        self.0.as_bstr()
    }
}

impl fmt::Debug for Array<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl<'a> crate::feature::record::attributes::field::value::Array<'a> for Array<'a> {
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Cow<'a, BStr>>> + 'a> {
        Box::new(self.iter().map(Ok))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iter() {
        let array = Array::new(b"8,13%2C21");

        assert_eq!(
            array.iter().collect::<Vec<_>>(),
            [Cow::from("8"), Cow::from("13,21")]
        );
    }
}
