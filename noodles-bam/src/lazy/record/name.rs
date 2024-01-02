use std::io;

use noodles_sam as sam;

/// A raw BAM record name.
#[derive(Debug, Eq, PartialEq)]
pub struct Name<'a>(&'a [u8]);

impl<'a> Name<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }

    /// Returns the name as a byte slice.
    ///
    /// The returned slice will _not_ have the trailing `NUL` terminator.
    pub fn as_bytes(&self) -> &[u8] {
        const NUL: u8 = 0x00;
        self.as_ref().strip_suffix(&[NUL]).unwrap_or(self.as_ref())
    }
}

impl<'a> sam::alignment::record::Name for Name<'a> {
    fn as_bytes(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl<'a> AsRef<[u8]> for Name<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryFrom<Name<'a>> for sam::record::Name {
    type Error = io::Error;

    fn try_from(bam_name: Name<'a>) -> Result<Self, Self::Error> {
        Self::try_from(bam_name.as_bytes().to_vec())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_as_bytes() {
        let name = Name::new(b"r0\x00");
        assert_eq!(name.as_bytes(), b"r0");

        let name = Name::new(b"r0");
        assert_eq!(name.as_bytes(), b"r0");
    }
}
