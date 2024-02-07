/// Raw VCF record filters.
#[derive(Debug, Eq, PartialEq)]
pub struct Filters<'a>(&'a str);

impl<'a> Filters<'a> {
    pub(super) fn new(buf: &'a str) -> Self {
        Self(buf)
    }

    /// Returns an iterator over all filters.
    pub fn iter(&self) -> impl Iterator<Item = &str> {
        const DELIMITER: char = ';';
        self.0.split(DELIMITER)
    }
}

impl<'a> AsRef<str> for Filters<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iter() {
        let filters = Filters::new("PASS");
        let actual: Vec<_> = filters.iter().collect();
        assert_eq!(actual, ["PASS"]);

        let filters = Filters::new("q10;s50");
        let actual: Vec<_> = filters.iter().collect();
        assert_eq!(actual, ["q10", "s50"]);
    }
}
