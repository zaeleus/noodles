use bstr::{BStr, ByteSlice};

/// A GFF directive.
pub struct Directive<'l> {
    src: &'l [u8],
    mid: usize,
}

impl<'l> Directive<'l> {
    pub(super) fn new(src: &'l [u8]) -> Self {
        let mid = src
            .iter()
            .position(|b| b.is_ascii_whitespace())
            .unwrap_or(src.len());

        Self { src, mid }
    }

    /// Returns the key.
    pub fn key(&self) -> &'l BStr {
        self.src[..self.mid].as_bstr()
    }

    /// Returns the value.
    pub fn value(&self) -> Option<&'l BStr> {
        self.src.get(self.mid + 1..).map(|s| s.as_bstr())
    }
}

impl AsRef<[u8]> for Directive<'_> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_key() {
        let directive = Directive::new(b"gff-version 3");
        assert_eq!(directive.key(), "gff-version");

        let directive = Directive::new(b"FASTA");
        assert_eq!(directive.key(), "FASTA");
    }

    #[test]
    fn test_value() {
        let directive = Directive::new(b"gff-version 3");
        assert_eq!(directive.value(), Some(BStr::new("3")));

        let directive = Directive::new(b"FASTA");
        assert!(directive.value().is_none());
    }
}
