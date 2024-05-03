//! Variant record samples.

pub mod keys;
pub mod sample;
mod series;

use std::io;

use self::sample::Value;
pub use self::{keys::Keys, sample::Sample, series::Series};
use crate::Header;

/// A variant record samples buffer.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Samples {
    pub(crate) keys: Keys,
    pub(crate) values: Vec<Vec<Option<Value>>>,
}

impl Samples {
    /// Creates a variant record samples buffer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::record_buf::{samples::Keys, Samples};
    /// let genotypes = Samples::new(Keys::default(), Vec::new());
    /// ```
    pub fn new(keys: Keys, values: Vec<Vec<Option<Value>>>) -> Self {
        Self { keys, values }
    }

    /// Returns whether there are any samples.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::record_buf::Samples;
    /// let samples = Samples::default();
    /// assert!(samples.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Returns the keys.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::{
    ///     record::samples::keys::key,
    ///     record_buf::{samples::Keys, Samples},
    /// };
    ///
    /// let samples = Samples::default();
    /// assert!(samples.keys().as_ref().is_empty());
    ///
    /// let keys: Keys = [String::from(key::GENOTYPE)].into_iter().collect();
    /// let samples = Samples::new(keys.clone(), Vec::new());
    /// assert_eq!(samples.keys(), &keys);
    /// ```
    pub fn keys(&self) -> &Keys {
        &self.keys
    }

    /// Returns a mutable reference to the keys.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::{
    ///     record::samples::keys::key,
    ///     record_buf::{samples::Keys, Samples},
    /// };
    ///
    /// let keys: Keys = [String::from(key::GENOTYPE)].into_iter().collect();
    ///
    /// let mut samples = Samples::default();
    /// *samples.keys_mut() = keys.clone();
    ///
    /// assert_eq!(samples.keys(), &keys);
    /// ```
    pub fn keys_mut(&mut self) -> &mut Keys {
        &mut self.keys
    }

    /// Returns samples.
    pub fn values(&self) -> impl Iterator<Item = Sample<'_>> {
        self.values
            .iter()
            .map(|values| Sample::new(&self.keys, values))
    }

    /// Returns the sample with the given sample name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     variant::{
    ///         record::samples::keys::key,
    ///         record_buf::{samples::sample::Value, Samples},
    ///     },
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_sample_name("sample0")
    ///     .add_sample_name("sample1")
    ///     .add_sample_name("sample2")
    ///     .build();
    ///
    /// let keys = [String::from(key::GENOTYPE)].into_iter().collect();
    /// let samples = Samples::new(
    ///     keys,
    ///     vec![
    ///         vec![Some(Value::from("0|0"))],
    ///         vec![Some(Value::from("1/1"))],
    ///         vec![],
    ///     ],
    /// );
    ///
    /// let sample = samples.get(&header, "sample0");
    /// assert_eq!(sample.and_then(|s| s.values().get(0)), Some(&Some(Value::from("0|0"))));
    /// ```
    pub fn get(&self, header: &Header, sample_name: &str) -> Option<Sample<'_>> {
        header
            .sample_names()
            .get_index_of(sample_name)
            .and_then(|i| self.get_index(i))
    }

    /// Returns the sample at the given index.
    pub fn get_index(&self, i: usize) -> Option<Sample<'_>> {
        self.values
            .get(i)
            .map(|values| Sample::new(&self.keys, values))
    }

    /// Returns the series with the given column name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::{
    ///     record::samples::keys::key,
    ///     record_buf::{samples::sample::Value, Samples},
    /// };
    ///
    /// let keys = [String::from(key::GENOTYPE)].into_iter().collect();
    /// let samples = Samples::new(
    ///     keys,
    ///     vec![
    ///         vec![Some(Value::from("0|0"))],
    ///         vec![Some(Value::from("1/1"))],
    ///         vec![],
    ///     ],
    /// );
    ///
    /// let series = samples.select(key::GENOTYPE).expect("missing genotype column");
    /// assert_eq!(series.get(0), Some(Some(&Value::from("0|0"))));
    /// assert_eq!(series.get(1), Some(Some(&Value::from("1/1"))));
    /// assert_eq!(series.get(2), Some(None));
    /// assert_eq!(series.get(3), None);
    /// ```
    pub fn select(&self, name: &str) -> Option<Series<'_>> {
        self.keys()
            .as_ref()
            .get_full(name)
            .map(|(i, name)| Series::new(name, &self.values[..], i))
    }

    /// Returns an iterator over series.
    pub fn series(&self) -> impl Iterator<Item = Series<'_>> {
        let column_names = self.keys();
        let column_count = column_names.as_ref().len();

        (0..column_count).map(|i| {
            // SAFETY: `i` is < `column_count`.
            let name = column_names.as_ref().get_index(i).unwrap();
            Series::new(name, &self.values[..], i)
        })
    }
}

impl crate::variant::record::Samples for Samples {
    fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    fn len(&self) -> usize {
        self.values.len()
    }

    fn column_names<'a, 'h: 'a>(
        &'a self,
        _: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<&'a str>> + 'a> {
        Box::new(self.keys.as_ref().iter().map(|key| Ok(key.as_str())))
    }

    fn select<'a, 'h: 'a>(
        &'a self,
        _: &'h Header,
        column_name: &str,
    ) -> Option<io::Result<Box<dyn crate::variant::record::samples::Series + 'a>>> {
        self.select(column_name)
            .map(|series| Box::new(series) as Box<dyn crate::variant::record::samples::Series>)
            .map(Ok)
    }

    fn series(
        &self,
    ) -> Box<
        dyn Iterator<Item = io::Result<Box<dyn crate::variant::record::samples::Series + '_>>> + '_,
    > {
        Box::new(
            self.series()
                .map(|series| Box::new(series) as Box<dyn crate::variant::record::samples::Series>)
                .map(Ok),
        )
    }

    fn iter(
        &self,
    ) -> Box<dyn Iterator<Item = Box<dyn crate::variant::record::samples::Sample + '_>> + '_> {
        Box::new(
            self.values()
                .map(|sample| Box::new(sample) as Box<dyn crate::variant::record::samples::Sample>),
        )
    }
}

impl crate::variant::record::Samples for &Samples {
    fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    fn len(&self) -> usize {
        self.values.len()
    }

    fn column_names<'a, 'h: 'a>(
        &'a self,
        _: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<&'a str>> + 'a> {
        Box::new(self.keys.as_ref().iter().map(|key| Ok(key.as_str())))
    }

    fn select<'a, 'h: 'a>(
        &'a self,
        _: &'h Header,
        column_name: &str,
    ) -> Option<io::Result<Box<dyn crate::variant::record::samples::Series + 'a>>> {
        Samples::select(self, column_name)
            .map(|series| Box::new(series) as Box<dyn crate::variant::record::samples::Series>)
            .map(Ok)
    }

    fn series(
        &self,
    ) -> Box<
        dyn Iterator<Item = io::Result<Box<dyn crate::variant::record::samples::Series + '_>>> + '_,
    > {
        Box::new(
            Samples::series(self)
                .map(|series| Box::new(series) as Box<dyn crate::variant::record::samples::Series>)
                .map(Ok),
        )
    }

    fn iter(
        &self,
    ) -> Box<dyn Iterator<Item = Box<dyn crate::variant::record::samples::Sample + '_>> + '_> {
        Box::new(
            self.values()
                .map(|sample| Box::new(sample) as Box<dyn crate::variant::record::samples::Sample>),
        )
    }
}

impl From<Samples> for (Keys, Vec<Vec<Option<Value>>>) {
    fn from(samples: Samples) -> Self {
        (samples.keys, samples.values)
    }
}
