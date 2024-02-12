mod keys;
mod sample;

use std::iter;

pub use self::{keys::Keys, sample::Sample};

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

    /// Returns the keys.
    pub fn keys(&self) -> Keys<'_> {
        let (src, _) = self.0.split_once(DELIMITER).unwrap_or_default();
        Keys::new(src)
    }

    /// Returns an iterator over samples.
    pub fn iter(&self) -> impl Iterator<Item = Option<Sample<'_>>> + '_ {
        let (_, mut src) = self.0.split_once(DELIMITER).unwrap_or_default();

        iter::from_fn(move || {
            if src.is_empty() {
                None
            } else {
                Some(parse_sample(&mut src, self.keys()))
            }
        })
    }
}

impl<'a> AsRef<str> for Samples<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

fn parse_sample<'a>(src: &mut &'a str, keys: Keys<'a>) -> Option<Sample<'a>> {
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
        _ => Some(Sample::new(buf, keys)),
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
    fn test_iter() {
        let samples = Samples::new("");
        assert!(samples.iter().next().is_none());

        let samples = Samples::new("GT:GQ\t0|0:13\t.");
        let actual: Vec<_> = samples.iter().collect();
        let expected = [Some(Sample::new("0|0:13", samples.keys())), None];
        assert_eq!(actual, expected);
    }
}
