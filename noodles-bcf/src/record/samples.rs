mod sample;
mod series;

use std::{io, iter};

use noodles_vcf as vcf;

use self::series::read_series;
pub use self::{sample::Sample, series::Series};
use crate::header::string_maps::StringStringMap;

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

    /// Converts BCF record genotypes to VCF record genotypes.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf::{header::string_maps::StringMap, record::Samples};
    /// use noodles_vcf as vcf;
    ///
    /// let bcf_genotypes = Samples::default();
    ///
    /// let header = vcf::Header::default();
    /// let string_maps = StringMap::default();
    /// let vcf_genotypes = bcf_genotypes.try_into_vcf_record_genotypes(&header, &string_maps)?;
    ///
    /// assert!(vcf_genotypes.is_empty());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn try_into_vcf_record_genotypes(
        &self,
        header: &vcf::Header,
        string_map: &StringStringMap,
    ) -> io::Result<vcf::record::Genotypes> {
        use crate::record::codec::decoder::read_genotypes;

        if self.is_empty() {
            return Ok(vcf::record::Genotypes::default());
        }

        let mut reader = self.src;

        let genotypes = read_genotypes(
            &mut reader,
            header.formats(),
            string_map,
            self.len(),
            self.format_count(),
        )
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        Ok(genotypes)
    }

    /// Returns the number of samples.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Samples;
    /// let samples = Samples::default();
    /// assert_eq!(samples.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.sample_count
    }

    /// Returns whether there are any samples.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Samples;
    /// let samples = Samples::default();
    /// assert!(samples.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.len() == 0
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
        (0..self.len()).map(|i| Sample::new(self, i))
    }
}

impl<'a> AsRef<[u8]> for Samples<'a> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}
