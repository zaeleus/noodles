use super::Keys;

/// A VCF record samples sample.
#[derive(Debug, Eq, PartialEq)]
pub struct Sample<'a> {
    src: &'a str,
    keys: Keys<'a>,
}

impl<'a> Sample<'a> {
    pub(super) fn new(src: &'a str, keys: Keys<'a>) -> Self {
        Self { src, keys }
    }

    pub fn values(&self) -> impl Iterator<Item = Option<&str>> + '_ {
        const MISSING: &str = ".";
        const DELIMITER: char = ':';

        self.src.split(DELIMITER).map(|s| match s {
            MISSING => None,
            _ => Some(s),
        })
    }

    /// Returns an iterator over fields.
    pub fn iter(&self) -> impl Iterator<Item = (&str, Option<&str>)> + '_ {
        self.keys.iter().zip(self.values())
    }
}

impl<'a> AsRef<str> for Sample<'a> {
    fn as_ref(&self) -> &str {
        self.src
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_values() {
        let keys = Keys::new("GT:GQ");
        let sample = Sample::new("0|0:.", keys);
        let actual: Vec<_> = sample.values().collect();
        let expected = [Some("0|0"), None];
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_iter() {
        let keys = Keys::new("GT:GQ");
        let sample = Sample::new("0|0:.", keys);
        let actual: Vec<_> = sample.iter().collect();
        let expected = [("GT", Some("0|0")), ("GQ", None)];
        assert_eq!(actual, expected);
    }
}
