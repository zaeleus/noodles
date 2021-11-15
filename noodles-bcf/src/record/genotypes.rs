use std::io;

use noodles_vcf as vcf;

use crate::header::StringMap;

/// BCF record genotypes.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Genotypes {
    buf: Vec<u8>,
    format_count: usize,
    sample_count: usize,
}

impl Genotypes {
    /// Converts BCF record genotypes to VCF record genotypes.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf::{header::StringMap, record::Genotypes};
    ///
    /// let bcf_genotypes = Genotypes::default();
    /// let string_map = StringMap::default();
    ///
    /// let (format, vcf_genotypes) = bcf_genotypes.try_into_vcf_record_genotypes(&string_map)?;
    ///
    /// assert!(format.is_none());
    /// assert!(vcf_genotypes.is_empty());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn try_into_vcf_record_genotypes(
        &self,
        string_map: &StringMap,
    ) -> io::Result<(Option<vcf::record::Format>, vcf::record::Genotypes)> {
        use crate::reader::record::read_genotypes;

        if self.is_empty() {
            return Ok((None, vcf::record::Genotypes::default()));
        }

        let mut reader = &self.buf[..];
        let genotypes = read_genotypes(&mut reader, string_map, self.len(), self.format_count())?;

        let first_genotype = genotypes.first().expect("unexpected empty genotypes");
        let keys: Vec<_> = first_genotype.keys().cloned().collect();
        let format = vcf::record::Format::try_from(keys)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        Ok((Some(format), genotypes))
    }

    /// Returns the number of samples.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Genotypes;
    /// let genotypes = Genotypes::default();
    /// assert_eq!(genotypes.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.sample_count
    }

    /// Returns whether there are any samples.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Genotypes;
    /// let genotypes = Genotypes::default();
    /// assert!(genotypes.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the number of fields per sample.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Genotypes;
    /// let genotypes = Genotypes::default();
    /// assert_eq!(genotypes.format_count(), 0);
    /// ```
    pub fn format_count(&self) -> usize {
        self.format_count
    }

    pub(crate) fn set_format_count(&mut self, format_count: usize) {
        self.format_count = format_count;
    }

    pub(crate) fn set_sample_count(&mut self, sample_count: usize) {
        self.sample_count = sample_count;
    }
}

impl AsRef<[u8]> for Genotypes {
    fn as_ref(&self) -> &[u8] {
        &self.buf
    }
}

impl AsMut<Vec<u8>> for Genotypes {
    fn as_mut(&mut self) -> &mut Vec<u8> {
        &mut self.buf
    }
}
