/// A raw VCF record genotypes sample.
#[derive(Debug, Eq, PartialEq)]
pub struct Sample<'a>(&'a str);

impl<'a> Sample<'a> {
    pub(super) fn new(buf: &'a str) -> Self {
        Self(buf)
    }

    pub fn values(&self) -> impl Iterator<Item = Option<&str>> + '_ {
        const MISSING: &str = ".";
        const DELIMITER: char = ':';

        self.0.split(DELIMITER).map(|s| match s {
            MISSING => None,
            _ => Some(s),
        })
    }
}

impl<'a> AsRef<str> for Sample<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_values() {
        let sample = Sample::new("0|0:.");
        let actual: Vec<_> = sample.values().collect();
        let expected = [Some("0|0"), None];
        assert_eq!(actual, expected);
    }
}
