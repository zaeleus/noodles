//! VCF record samples series.

pub mod value;

use std::io;

use super::Samples;
use crate::{Header, variant::record::samples::series::Value};

/// VCF record samples series.
pub struct Series<'r> {
    name: &'r str,
    samples: &'r Samples<'r>,
    i: usize,
}

impl<'r> Series<'r> {
    pub(super) fn new(name: &'r str, samples: &'r Samples<'r>, i: usize) -> Self {
        Self { name, samples, i }
    }
}

impl crate::variant::record::samples::Series for Series<'_> {
    fn name<'a, 'h: 'a>(&'a self, _: &'h Header) -> io::Result<&'a str> {
        Ok(self.name)
    }

    fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        i: usize,
    ) -> Option<Option<io::Result<Value<'a>>>> {
        let sample = self.samples.iter().nth(i)?;
        sample.get_index(header, self.i)
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<Option<Value<'a>>>> + 'a> {
        Box::new(
            self.samples
                .iter()
                .map(|sample| match sample.get_index(header, self.i) {
                    Some(value) => value.transpose(),
                    None => Ok(None),
                }),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record::{samples::Series as _, samples::keys::key};

    #[test]
    fn test_name() {
        let header = Header::default();

        let samples = Samples::new("GT:GQ\t0|0:13\t0/1:8");
        let series = Series::new(key::CONDITIONAL_GENOTYPE_QUALITY, &samples, 1);

        assert!(matches!(
            series.name(&header),
            Ok(name) if name == key::CONDITIONAL_GENOTYPE_QUALITY
        ));
    }

    #[test]
    fn test_get() {
        let header = Header::default();

        let samples = Samples::new("GT:GQ\t0|0:13\t0/1:8");
        let series = Series::new(key::CONDITIONAL_GENOTYPE_QUALITY, &samples, 1);

        assert!(matches!(
            series.get(&header, 0),
            Some(Some(Ok(Value::Integer(13))))
        ));

        assert!(matches!(
            series.get(&header, 1),
            Some(Some(Ok(Value::Integer(8))))
        ));

        assert!(series.get(&header, 2).is_none());
    }
}
