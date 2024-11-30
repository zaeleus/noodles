/// A GFF directive.
pub struct Directive<'l> {
    src: &'l str,
    mid: usize,
}

impl<'l> Directive<'l> {
    pub(super) fn new(src: &'l str) -> Self {
        let mid = src
            .find(|c: char| c.is_ascii_whitespace())
            .unwrap_or(src.len());

        Self { src, mid }
    }

    /// Returns the key.
    pub fn key(&self) -> &'l str {
        &self.src[..self.mid]
    }

    /// Returns the value.
    pub fn value(&self) -> Option<&'l str> {
        self.src.get(self.mid + 1..)
    }
}

impl AsRef<str> for Directive<'_> {
    fn as_ref(&self) -> &str {
        self.src
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_key() {
        let directive = Directive::new("gff-version 3");
        assert_eq!(directive.key(), "gff-version");

        let directive = Directive::new("FASTA");
        assert_eq!(directive.key(), "FASTA");
    }

    #[test]
    fn test_value() {
        let directive = Directive::new("gff-version 3");
        assert_eq!(directive.value(), Some("3"));

        let directive = Directive::new("FASTA");
        assert!(directive.value().is_none());
    }
}
