use std::io;

use noodles_vcf as vcf;

use crate::header::string_maps::StringStringMap;

/// BCF record filters.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Filters(Vec<usize>);

impl Filters {
    /// Converts BCF record filters to VCF record filters.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf::{header::StringMaps, lazy::record::Filters};
    ///
    /// let bcf_filters = Filters::default();
    /// let string_maps = StringMaps::default();
    /// let vcf_filters = bcf_filters.try_into_vcf_record_filters(string_maps.strings())?;
    ///
    /// assert!(vcf_filters.is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn try_into_vcf_record_filters(
        &self,
        string_string_map: &StringStringMap,
    ) -> io::Result<Option<vcf::record::Filters>> {
        let raw_filters: Vec<_> = self
            .0
            .iter()
            .map(|&i| {
                string_string_map.get_index(i).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid string map index: {i}"),
                    )
                })
            })
            .collect::<Result<_, _>>()?;

        if raw_filters.is_empty() {
            Ok(None)
        } else {
            vcf::record::Filters::try_from_iter(raw_filters)
                .map(Some)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        }
    }

    /// Returns the number of filters.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::lazy::record::Filters;
    /// let filters = Filters::default();
    /// assert_eq!(filters.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns whether there are any filters.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::lazy::record::Filters;
    /// let filters = Filters::default();
    /// assert!(filters.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Removes all filter IDs from the filters list.
    ///
    /// This does not affect the capacity of the list.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::lazy::record::Filters;
    /// let mut filters = Filters::default();
    /// filters.clear();
    /// assert!(filters.is_empty());
    /// ```
    pub fn clear(&mut self) {
        self.0.clear();
    }
}

impl AsRef<[usize]> for Filters {
    fn as_ref(&self) -> &[usize] {
        &self.0
    }
}

impl AsMut<Vec<usize>> for Filters {
    fn as_mut(&mut self) -> &mut Vec<usize> {
        &mut self.0
    }
}
