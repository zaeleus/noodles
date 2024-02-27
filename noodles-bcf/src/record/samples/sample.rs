use std::io;

use noodles_vcf::{self as vcf, variant::record::samples::series::Value};

use super::Samples;

/// A BCF record sample.
pub struct Sample<'r> {
    samples: &'r Samples<'r>,
    i: usize,
}

impl<'r> Sample<'r> {
    pub(super) fn new(samples: &'r Samples<'r>, i: usize) -> Self {
        Self { samples, i }
    }
}

impl<'r> vcf::variant::record::samples::Sample for Sample<'r> {
    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
    ) -> Box<dyn Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a> {
        Box::new(self.samples.series().map(|result| {
            result.and_then(|series| {
                let name = series.name(header)?;

                let value = series
                    .get(self.i)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing value"))?;

                Ok((name, value))
            })
        }))
    }
}
