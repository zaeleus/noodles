use std::io;

/// Raw GFF record strand.
pub struct Strand<'a>(&'a str);

impl<'a> Strand<'a> {
    pub(super) fn new(buf: &'a str) -> Self {
        Self(buf)
    }
}

impl<'a> AsRef<str> for Strand<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

impl<'a> TryFrom<Strand<'a>> for crate::record::Strand {
    type Error = io::Error;

    fn try_from(raw_strand: Strand<'a>) -> Result<Self, Self::Error> {
        raw_strand
            .as_ref()
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}
