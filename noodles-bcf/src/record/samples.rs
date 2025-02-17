//! BCF record samples.

mod sample;
pub mod series;

use std::{io, iter};

use noodles_vcf as vcf;

use self::series::read_series;
pub use self::{sample::Sample, series::Series};

/// BCF record samples.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Samples<'r> {
    src: &'r [u8],
    sample_count: usize,
    format_count: usize,
}

impl<'r> Samples<'r> {
    pub(super) fn new(src: &'r [u8], sample_count: usize, format_count: usize) -> Self {
        Self {
            src,
            sample_count,
            format_count,
        }
    }

    /// Returns the number of fields per sample.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Samples;
    /// let samples = Samples::default();
    /// assert_eq!(samples.format_count(), 0);
    /// ```
    pub fn format_count(&self) -> usize {
        self.format_count
    }

    /// Returns the sample with the given sample name.
    pub fn get(&'r self, header: &vcf::Header, sample_name: &str) -> Option<Sample<'r>> {
        header
            .sample_names()
            .get_index_of(sample_name)
            .and_then(|i| self.get_index(i))
    }

    /// Returns a sample at the given index.
    pub fn get_index(&'r self, i: usize) -> Option<Sample<'r>> {
        if i < self.sample_count {
            Some(Sample::new(self, i))
        } else {
            None
        }
    }

    /// Returns the series with the given column name.
    pub fn select<'h: 'r>(
        &'r self,
        header: &'h vcf::Header,
        column_name: &str,
    ) -> Option<io::Result<Series<'r>>> {
        for result in self.series() {
            let series = match result {
                Ok(s) => s,
                Err(e) => return Some(Err(e)),
            };

            match series.name(header) {
                Ok(name) if name == column_name => return Some(Ok(series)),
                Ok(_) => {}
                Err(e) => return Some(Err(e)),
            }
        }

        None
    }

    /// Returns an iterator over series.
    pub fn series(&'r self) -> impl Iterator<Item = io::Result<Series<'r>>> + 'r {
        let mut src = self.src;

        iter::from_fn(move || {
            if src.is_empty() {
                None
            } else {
                Some(read_series(&mut src, self.sample_count))
            }
        })
    }

    /// Returns an iterator over samples.
    pub fn iter(&self) -> impl Iterator<Item = Sample<'_>> {
        (0..self.sample_count).map(|i| Sample::new(self, i))
    }
}

impl AsRef<[u8]> for Samples<'_> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}

impl vcf::variant::record::Samples for Samples<'_> {
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn len(&self) -> usize {
        self.sample_count
    }

    fn column_names<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
    ) -> Box<dyn Iterator<Item = io::Result<&'a str>> + 'a> {
        Box::new(
            self.series()
                .map(|result| result.and_then(|series| series.name(header))),
        )
    }

    fn select<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
        column_name: &str,
    ) -> Option<io::Result<Box<dyn vcf::variant::record::samples::Series + 'a>>> {
        self.select(header, column_name).map(|result| {
            result.map(|series| Box::new(series) as Box<dyn vcf::variant::record::samples::Series>)
        })
    }

    fn series(
        &self,
    ) -> Box<
        dyn Iterator<Item = io::Result<Box<dyn vcf::variant::record::samples::Series + '_>>> + '_,
    > {
        Box::new(self.series().map(|result| {
            result.map(|series| Box::new(series) as Box<dyn vcf::variant::record::samples::Series>)
        }))
    }

    fn iter(
        &self,
    ) -> Box<dyn Iterator<Item = Box<dyn vcf::variant::record::samples::Sample + '_>> + '_> {
        Box::new(
            self.iter()
                .map(|sample| Box::new(sample) as Box<dyn vcf::variant::record::samples::Sample>),
        )
    }
}
