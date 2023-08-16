mod sample;

use std::{io, iter};

pub use self::sample::Sample;
use crate::record::{FIELD_DELIMITER, MISSING_FIELD};

/// Raw VCF record genotypes.
#[derive(Debug, Eq, PartialEq)]
pub struct Genotypes<'a>(&'a str);

impl<'a> Genotypes<'a> {
    pub(super) fn new(buf: &'a str) -> Self {
        Self(buf)
    }

    /// Returns whether there may be any genotypes.
    pub fn is_empty(&self) -> bool {
        let is_missing = self
            .0
            .split(FIELD_DELIMITER)
            .next()
            .map(|s| s == MISSING_FIELD)
            .unwrap_or_default();

        self.0.is_empty() || is_missing
    }

    /// Returns an iterator over keys.
    pub fn keys(&self) -> io::Result<Box<dyn Iterator<Item = &str> + '_>> {
        const DELIMITER: char = ':';

        if self.is_empty() {
            return Ok(Box::new(iter::empty()));
        }

        let (raw_format, _) = self
            .0
            .split_once(FIELD_DELIMITER)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing field separator"))?;

        Ok(Box::new(raw_format.split(DELIMITER)))
    }

    /// Returns an iterator over samples.
    pub fn samples(&self) -> io::Result<Box<dyn Iterator<Item = Option<Sample<'_>>> + '_>> {
        if self.is_empty() {
            return Ok(Box::new(iter::empty()));
        }

        let (_, raw_samples) = self
            .0
            .split_once(FIELD_DELIMITER)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing field separator"))?;

        Ok(Box::new(raw_samples.split(FIELD_DELIMITER).map(
            |s| match s {
                "." => None,
                _ => Some(Sample::new(s)),
            },
        )))
    }
}

impl<'a> AsRef<str> for Genotypes<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        assert!(Genotypes::new("").is_empty());
        assert!(Genotypes::new(".\t.").is_empty());
        assert!(!Genotypes::new("GT:GQ\t0|0:13").is_empty());
    }

    #[test]
    fn test_keys() -> io::Result<()> {
        let genotypes = Genotypes::new("");
        assert!(genotypes.keys()?.next().is_none());

        let genotypes = Genotypes::new(".\t.");
        assert!(genotypes.keys()?.next().is_none());

        let genotypes = Genotypes::new("GT:GQ\t0|0:13");
        let actual: Vec<_> = genotypes.keys()?.collect();
        let expected = ["GT", "GQ"];
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_samples() -> io::Result<()> {
        let genotypes = Genotypes::new("");
        assert!(genotypes.samples()?.next().is_none());

        let genotypes = Genotypes::new(".\t.");
        assert!(genotypes.samples()?.next().is_none());

        let genotypes = Genotypes::new("GT:GQ\t0|0:13\t.");
        let actual: Vec<_> = genotypes.samples()?.collect();
        let expected = [Some(Sample::new("0|0:13")), None];
        assert_eq!(actual, expected);

        Ok(())
    }
}
