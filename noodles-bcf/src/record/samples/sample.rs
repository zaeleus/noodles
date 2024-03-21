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
    fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
        key: &str,
    ) -> Option<io::Result<Option<Value<'a>>>> {
        for result in self.iter(header) {
            match result {
                Ok((k, v)) => {
                    if k == key {
                        return Some(Ok(v));
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        None
    }

    fn get_index<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
        i: usize,
    ) -> Option<io::Result<Option<Value<'a>>>> {
        self.iter(header)
            .nth(i)
            .map(|result| result.map(|(_, value)| value))
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
    ) -> Box<dyn Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a> {
        let series = self.samples.series();

        Box::new(series.map(|result| {
            result.and_then(|series| {
                let name = series.name(header)?;

                let value = series
                    .get(header, self.i)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing value"))?
                    .transpose()?;

                Ok((name, value))
            })
        }))
    }
}
