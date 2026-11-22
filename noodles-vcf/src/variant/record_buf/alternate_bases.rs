use std::io;

/// A variant record alternate bases buffer.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AlternateBases(Vec<String>);

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

impl Extend<String> for AlternateBases {
    fn extend<T: IntoIterator<Item = String>>(&mut self, iter: T) {
        self.0.extend(iter);
    }
}

impl From<Vec<String>> for AlternateBases {
    fn from(alleles: Vec<String>) -> Self {
        Self(alleles)
    }
}

impl FromIterator<String> for AlternateBases {
    fn from_iter<T: IntoIterator<Item = String>>(iter: T) -> Self {
        let mut alternate_bases = Self::default();
        alternate_bases.extend(iter);
        alternate_bases
    }
}

impl crate::variant::record::AlternateBases for AlternateBases {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<&str>>> + '_> {
        Box::new(self.0.iter().map(|allele| Some(allele.as_ref())).map(Ok))
    }
}

impl crate::variant::record::AlternateBases for &AlternateBases {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<&str>>> + '_> {
        Box::new(self.0.iter().map(|allele| Some(allele.as_ref())).map(Ok))
    }
}
