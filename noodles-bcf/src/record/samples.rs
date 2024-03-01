mod sample;
mod series;

use std::{io, iter};

use noodles_vcf::{self as vcf, variant::record::Samples as _};

use self::series::read_series;
pub use self::{sample::Sample, series::Series};

/// BCF record genotypes.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Samples<'a> {
    src: &'a [u8],
    sample_count: usize,
    format_count: usize,
}

impl<'a> Samples<'a> {
    pub(super) fn new(src: &'a [u8], sample_count: usize, format_count: usize) -> Self {
        Self {
            src,
            sample_count,
            format_count,
        }
    }

    /// Converts BCF record samples to VCF record samples.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf::record::Samples;
    /// use noodles_vcf as vcf;
    ///
    /// let bcf_samples = Samples::default();
    ///
    /// let header = vcf::Header::default();
    /// let vcf_samples = bcf_samples.try_into_vcf_record_samples(&header)?;
    ///
    /// assert!(vcf_samples.is_empty());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn try_into_vcf_record_samples(
        &self,
        header: &vcf::Header,
    ) -> io::Result<vcf::variant::record_buf::Samples> {
        use crate::record::codec::decoder::read_samples;

        if self.is_empty() {
            return Ok(vcf::variant::record_buf::Samples::default());
        }

        let mut reader = self.src;

        let genotypes = read_samples(&mut reader, header, self.len(), self.format_count())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        Ok(genotypes)
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

    /// Returns a sample at the given index.
    pub fn get_index(&self, i: usize) -> Option<Sample<'_>> {
        if i < self.sample_count {
            Some(Sample::new(self, i))
        } else {
            None
        }
    }

    /// Returns an iterator over series.
    pub fn series(&self) -> impl Iterator<Item = io::Result<Series<'_>>> {
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

impl<'a> AsRef<[u8]> for Samples<'a> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}

impl<'a> vcf::variant::record::Samples for Samples<'a> {
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn len(&self) -> usize {
        self.sample_count
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
