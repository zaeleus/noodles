use noodles_sam as sam;

/// A BAM record name.
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

impl<'a> From<Name<'a>> for sam::alignment::record_buf::Name {
    fn from(name: Name<'a>) -> Self {
        Self::from(name.as_ref())
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
