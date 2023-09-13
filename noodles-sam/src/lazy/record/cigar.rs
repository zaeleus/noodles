use std::io;

/// Raw SAM record CIGAR operations.
#[derive(Debug, Eq, PartialEq)]
pub struct Cigar<'a>(&'a [u8]);

impl<'a> Cigar<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }

    /// Returns whether there are any CIGAR operations.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl<'a> AsRef<[u8]> for Cigar<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryFrom<Cigar<'a>> for crate::record::Cigar {
    type Error = io::Error;

    fn try_from(Cigar(src): Cigar<'a>) -> Result<Self, Self::Error> {
        use crate::reader::record::parse_cigar;

        let mut cigar = crate::record::Cigar::default();

        if !src.is_empty() {
            parse_cigar(src, &mut cigar)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        }

        Ok(cigar)
    }
}
