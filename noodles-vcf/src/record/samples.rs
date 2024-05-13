//! VCF record samples.

mod keys;
mod sample;
pub mod series;

use std::{io, iter};

pub use self::{keys::Keys, sample::Sample, series::Series};
use crate::Header;

const DELIMITER: char = '\t';

/// Raw VCF record genotypes.
#[derive(Debug, Eq, PartialEq)]
pub struct Samples<'r>(&'r str);

impl<'r> Samples<'r> {
    pub(super) fn new(buf: &'r str) -> Self {
        Self(buf)
    }

    /// Returns whether there may be any genotypes.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the keys.
    pub fn keys(&self) -> Keys<'r> {
        let (src, _) = self.0.split_once(DELIMITER).unwrap_or_default();
        Keys::new(src)
    }

    /// Returns the sample with the given sample name.
    pub fn get(&self, header: &Header, sample_name: &str) -> Option<Sample<'r>> {
        header
            .sample_names()
            .get_index_of(sample_name)
            .and_then(|i| self.get_index(i))
    }

    /// Returns the sample at the given index.
    pub fn get_index(&self, i: usize) -> Option<Sample<'r>> {
        self.iter().nth(i)
    }

    /// Returns the series with the given column name.
    pub fn select(&'r self, column_name: &str) -> Option<Series<'r>> {
        self.keys()
            .iter()
            .enumerate()
            .find(|(_, key)| *key == column_name)
            .map(|(i, key)| Series::new(key, self, i))
    }

    /// Returns an iterator over series.
    pub fn series(&'r self) -> impl Iterator<Item = Series<'r>> + '_ {
        self.keys()
            .iter()
            .enumerate()
            .map(|(i, key)| Series::new(key, self, i))
    }

    /// Returns an iterator over samples.
    pub fn iter(&self) -> impl Iterator<Item = Sample<'r>> + '_ {
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

impl<'r> crate::variant::record::Samples for Samples<'r> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.iter().count()
    }

    fn column_names<'a, 'h: 'a>(
        &'a self,
        _: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<&'a str>> + 'a> {
        Box::new(self.keys().iter().map(Ok))
    }

    fn select<'a, 'h: 'a>(
        &'a self,
        _: &'h Header,
        column_name: &str,
    ) -> Option<io::Result<Box<dyn crate::variant::record::samples::Series + 'a>>> {
        self.select(column_name)
            .map(|series| Box::new(series) as Box<dyn crate::variant::record::samples::Series>)
            .map(Ok)
    }

    fn series(
        &self,
    ) -> Box<
        dyn Iterator<Item = io::Result<Box<dyn crate::variant::record::samples::Series + '_>>> + '_,
    > {
        Box::new(
            self.series()
                .map(|series| Box::new(series) as Box<dyn crate::variant::record::samples::Series>)
                .map(Ok),
        )
    }

    fn iter(
        &self,
    ) -> Box<dyn Iterator<Item = Box<dyn crate::variant::record::samples::Sample + '_>> + '_> {
        Box::new(
            self.iter()
                .map(|sample| Box::new(sample) as Box<dyn crate::variant::record::samples::Sample>),
        )
    }
}

fn parse_sample<'r>(src: &mut &'r str, keys: Keys<'r>) -> Sample<'r> {
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

    if buf == MISSING {
        Sample::new("", keys)
    } else {
        Sample::new(buf, keys)
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
    fn test_get() {
        let header = Header::builder()
            .add_sample_name("sample0")
            .add_sample_name("sample1")
            .add_sample_name("sample2")
            .build();

        let samples = Samples::new("GT\t0|0\t1/1\t.");

        let actual = samples.get(&header, "sample0");
        let expected = Sample::new("0|0", samples.keys());
        assert_eq!(actual, Some(expected));

        assert!(samples.get(&header, "sample3").is_none());
    }

    #[test]
    fn test_get_index() {
        let samples = Samples::new("GT\t0|0\t1/1\t.");
        let actual = samples.get_index(0);
        let expected = Sample::new("0|0", samples.keys());
        assert_eq!(actual, Some(expected));

        assert!(samples.get_index(3).is_none());
    }

    #[test]
    fn test_select() {
        use crate::variant::record::samples::keys::key;

        let samples = Samples::new("");
        assert!(samples.select(key::CONDITIONAL_GENOTYPE_QUALITY).is_none());

        let samples = Samples::new("GT:GQ\t0|0:13\t.");
        assert!(samples.select(key::CONDITIONAL_GENOTYPE_QUALITY).is_some());
    }

    #[test]
    fn test_iter() {
        let samples = Samples::new("");
        assert!(samples.iter().next().is_none());

        let samples = Samples::new("GT:GQ\t0|0:13\t.");
        let actual: Vec<_> = samples.iter().collect();
        let expected = [
            Sample::new("0|0:13", samples.keys()),
            Sample::new("", samples.keys()),
        ];
        assert_eq!(actual, expected);
    }
}
