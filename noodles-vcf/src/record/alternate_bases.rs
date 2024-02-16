use std::fmt;

use crate::variant::record::AlternateBases as _;

/// VCF record alternate bases.
pub struct AlternateBases<'a>(&'a str);

const DELIMITER: char = ',';

impl<'a> AlternateBases<'a> {
    pub(super) fn new(src: &'a str) -> Self {
        Self(src)
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

impl<'a> crate::variant::record::AlternateBases for AlternateBases<'a> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        if self.is_empty() {
            0
        } else {
            self.iter().count()
        }
    }

    fn iter(&self) -> Box<dyn Iterator<Item = &str> + '_> {
        Box::new(self.0.split(DELIMITER))
    }
}
