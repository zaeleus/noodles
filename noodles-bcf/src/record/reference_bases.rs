use std::io;

use noodles_vcf as vcf;

/// BCF record reference bases.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceBases<'a>(&'a [u8]);

impl<'a> ReferenceBases<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }
}

impl<'a> AsRef<[u8]> for ReferenceBases<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> vcf::variant::record::ReferenceBases for ReferenceBases<'a> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u8>> + '_> {
        Box::new(self.0.iter().copied().map(Ok))
    }
}
