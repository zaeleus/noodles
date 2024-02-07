mod sample;

use std::iter;

pub use self::sample::Sample;

const DELIMITER: char = '\t';

/// Raw VCF record genotypes.
#[derive(Debug, Eq, PartialEq)]
pub struct Samples<'a>(&'a str);

impl<'a> Samples<'a> {
    pub(super) fn new(buf: &'a str) -> Self {
        Self(buf)
    }

    /// Returns whether there may be any genotypes.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns an iterator over keys.
    pub fn keys(&self) -> impl Iterator<Item = &str> + '_ {
        let (mut src, _) = self.0.split_once(DELIMITER).unwrap_or_default();

        iter::from_fn(move || {
            if src.is_empty() {
                None
            } else {
                Some(parse_key(&mut src))
            }
        })
    }

    /// Returns an iterator over samples.
    pub fn samples(&self) -> impl Iterator<Item = Option<Sample<'_>>> + '_ {
        let (_, mut src) = self.0.split_once(DELIMITER).unwrap_or_default();

        iter::from_fn(move || {
            if src.is_empty() {
                None
            } else {
                Some(parse_sample(&mut src))
            }
        })
    }
}

impl<'a> AsRef<str> for Samples<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

fn parse_key<'a>(src: &mut &'a str) -> &'a str {
    const DELIMITER: u8 = b':';

    match src.as_bytes().iter().position(|&b| b == DELIMITER) {
        Some(i) => {
            let (buf, rest) = src.split_at(i);
            *src = &rest[1..];
            buf
        }
        None => {
            let (buf, rest) = src.split_at(src.len());
            *src = rest;
            buf
        }
    }
}

fn parse_sample<'a>(src: &mut &'a str) -> Option<Sample<'a>> {
    const DELIMITER: u8 = b'\t';
    const MISSING: &str = ".";

    let buf = match src.as_bytes().iter().position(|&b| b == DELIMITER) {
        Some(i) => {
            let (buf, rest) = src.split_at(i);
            *src = &rest[1..];
            buf
        }
        None => {
            let (buf, rest) = src.split_at(src.len());
            *src = rest;
            buf
        }
    };

    match buf {
        MISSING => None,
        _ => Some(Sample::new(buf)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        assert!(Samples::new("").is_empty());
        assert!(!Samples::new("GT:GQ\t0|0:13").is_empty());
    }

    #[test]
    fn test_keys() {
        let genotypes = Samples::new("");
        assert!(genotypes.keys().next().is_none());

        let genotypes = Samples::new("GT:GQ\t0|0:13");
        let actual: Vec<_> = genotypes.keys().collect();
        let expected = ["GT", "GQ"];
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_samples() {
        let genotypes = Samples::new("");
        assert!(genotypes.samples().next().is_none());

        let genotypes = Samples::new("GT:GQ\t0|0:13\t.");
        let actual: Vec<_> = genotypes.samples().collect();
        let expected = [Some(Sample::new("0|0:13")), None];
        assert_eq!(actual, expected);
    }
}
