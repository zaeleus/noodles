/// GTF record attributes.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Attributes(Vec<(String, String)>);

impl Attributes {
    /// Returns whether there are any entries.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf::record_buf::Attributes;
    /// let attributes = Attributes::default();
    /// assert!(attributes.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of entries.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf::record_buf::Attributes;
    /// let attributes = Attributes::default();
    /// assert_eq!(attributes.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns the value at the given key.
    ///
    /// This returns the first match of the key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf::record_buf::Attributes;
    ///
    /// let attributes: Attributes = [(String::from("id"), String::from("g0"))]
    ///     .into_iter()
    ///     .collect();
    ///
    /// assert_eq!(attributes.get("id"), Some("g0"));
    /// assert!(attributes.get("source").is_none());
    /// ```
    pub fn get<'a>(&'a self, key: &str) -> Option<&'a str> {
        self.0
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }
}

impl AsRef<[(String, String)]> for Attributes {
    fn as_ref(&self) -> &[(String, String)] {
        &self.0
    }
}

impl From<Vec<(String, String)>> for Attributes {
    fn from(entries: Vec<(String, String)>) -> Self {
        Self(entries)
    }
}

impl Extend<(String, String)> for Attributes {
    fn extend<T: IntoIterator<Item = (String, String)>>(&mut self, iter: T) {
        self.0.extend(iter);
    }
}

impl FromIterator<(String, String)> for Attributes {
    fn from_iter<T: IntoIterator<Item = (String, String)>>(iter: T) -> Self {
        let mut attributes = Self::default();
        attributes.extend(iter);
        attributes
    }
}
