//! Inner VCF header contig map value.

mod builder;
pub mod name;
pub(crate) mod tag;

pub use self::{name::Name, tag::Tag};

use super::{Indexed, Inner, Map};

/// An inner VCF header contig map value.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Contig {
    pub(crate) length: Option<usize>,
    pub(crate) md5: Option<String>,
    pub(crate) url: Option<String>,
    pub(crate) idx: Option<usize>,
}

impl Inner for Contig {
    type StandardTag = tag::Standard;
    type Builder = builder::Builder;
}

impl Indexed for Contig {
    fn idx(&self) -> Option<usize> {
        self.idx
    }

    fn idx_mut(&mut self) -> &mut Option<usize> {
        &mut self.idx
    }
}

impl Map<Contig> {
    /// Creates a VCF header contig map value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    /// let map = Map::<Contig>::new();
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    /// let map = Map::<Contig>::new();
    /// assert!(map.length().is_none());
    /// ```
    pub fn length(&self) -> Option<usize> {
        self.inner.length
    }

    /// Returns a mutable reference to the length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    ///
    /// let mut map = Map::<Contig>::new();
    /// assert!(map.length().is_none());
    ///
    /// *map.length_mut() = Some(8);
    /// assert_eq!(map.length(), Some(8));
    /// ```
    pub fn length_mut(&mut self) -> &mut Option<usize> {
        &mut self.inner.length
    }

    /// Returns the MD5 hexdigest.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    /// let map = Map::<Contig>::new();
    /// assert!(map.md5().is_none());
    /// ```
    pub fn md5(&self) -> Option<&str> {
        self.inner.md5.as_deref()
    }

    /// Returns a mutable reference to the MD5 hexdigest.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    ///
    /// let mut map = Map::<Contig>::new();
    /// assert!(map.md5().is_none());
    ///
    /// *map.md5_mut() = Some(String::from("d7eba311421bbc9d3ada44709dd61534"));
    /// assert_eq!(map.md5(), Some("d7eba311421bbc9d3ada44709dd61534"));
    /// ```
    pub fn md5_mut(&mut self) -> &mut Option<String> {
        &mut self.inner.md5
    }

    /// Returns the URL.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    /// let map = Map::<Contig>::new();
    /// assert!(map.url().is_none());
    /// ```
    pub fn url(&self) -> Option<&str> {
        self.inner.url.as_deref()
    }

    /// Returns a mutable reference to the URL.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::value::{map::Contig, Map};
    ///
    /// let mut map = Map::<Contig>::new();
    /// assert!(map.url().is_none());
    ///
    /// *map.url_mut() = Some(String::from("https://example.com/reference.fa"));
    /// assert_eq!(map.url(), Some("https://example.com/reference.fa"));
    /// ```
    pub fn url_mut(&mut self) -> &mut Option<String> {
        &mut self.inner.url
    }
}
