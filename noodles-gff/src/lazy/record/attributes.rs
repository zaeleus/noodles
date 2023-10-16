/// Raw GFF record attributes.
pub struct Attributes<'a>(&'a str);

impl<'a> Attributes<'a> {
    pub(super) fn new(buf: &'a str) -> Self {
        Self(buf)
    }

    /// Returns whether there are any attributes.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl<'a> AsRef<str> for Attributes<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        let attributes = Attributes::new("");
        assert!(attributes.is_empty());

        let attributes = Attributes::new("gene_id=ndls0;gene_name=gene0");
        assert!(!attributes.is_empty());
    }
}
