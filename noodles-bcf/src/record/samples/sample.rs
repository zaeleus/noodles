use std::io;

use noodles_vcf::{header::StringMaps, variant::record::samples::series::Value};

use super::Samples;

/// A BCF record sample.
pub struct Sample<'a> {
    samples: &'a Samples<'a>,
    i: usize,
}

impl<'a> Sample<'a> {
    pub(super) fn new(samples: &'a Samples<'a>, i: usize) -> Self {
        Self { samples, i }
    }

    /// Returns an iterator over fields.
    pub fn iter<'h>(
        &self,
        string_maps: &'h StringMaps,
    ) -> impl Iterator<Item = io::Result<(&'h str, Option<Value<'_>>)>> {
        self.samples.series().map(|result| {
            result.and_then(|series| {
                let name = series.name(string_maps)?;

                let value = series
                    .get(self.i)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing value"))?;

                Ok((name, value))
            })
        })
    }
}
