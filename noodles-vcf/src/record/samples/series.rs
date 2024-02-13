use std::io;

use super::Samples;
use crate::{variant::record::samples::series::Value, Header};

pub struct Series<'a> {
    name: &'a str,
    samples: &'a Samples<'a>,
    i: usize,
}

impl<'a> Series<'a> {
    pub(super) fn new(name: &'a str, samples: &'a Samples<'a>, i: usize) -> Self {
        Self { name, samples, i }
    }

    /// Returns the name.
    pub fn name(&self) -> &str {
        self.name
    }

    /// Returns the value at the given index.
    pub fn get<'h: 'a>(
        &self,
        header: &'h Header,
        i: usize,
    ) -> Option<Option<io::Result<Value<'a>>>> {
        let sample = self.samples.iter().nth(i)??;
        sample.get_index(header, self.i)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record_buf::samples::keys::key;

    #[test]
    fn test_name() {
        let samples = Samples::new("GT:GQ\t0|0:13\t0/1:8");
        let series = Series::new(key::CONDITIONAL_GENOTYPE_QUALITY, &samples, 1);
        assert_eq!(series.name(), key::CONDITIONAL_GENOTYPE_QUALITY);
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
