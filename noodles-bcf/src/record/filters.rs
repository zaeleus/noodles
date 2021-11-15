/// BCF record filters.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Filters(Vec<usize>);

impl Filters {
    /// Returns the number of filters.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Filters;
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
    /// use noodles_bcf::record::Filters;
    /// let filters = Filters::default();
    /// assert!(filters.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.len() == 0
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
