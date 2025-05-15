use std::collections::HashMap;

/// An indexed map of VCF strings.
///
/// This is also called a dictionary of strings.
///
/// See § 6.2.1 Dictionary of strings (2021-05-13).
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct StringMap {
    pub(super) indices: HashMap<String, usize>,
    pub(super) entries: Vec<Option<String>>,
}

impl StringMap {
    /// Returns an entry by index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::string_maps::StringMap;
    /// let string_map = StringMap::default();
    /// assert!(string_map.get_index(0).is_none());
    /// ```
    pub fn get_index(&self, i: usize) -> Option<&str> {
        self.entries.get(i).and_then(|entry| entry.as_deref())
    }

    /// Returns the index of the entry of the given value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::string_maps::StringMap;
    /// let string_map = StringMap::default();
    /// assert!(string_map.get_index_of("PASS").is_none());
    /// ```
    pub fn get_index_of(&self, value: &str) -> Option<usize> {
        self.indices.get(value).copied()
    }

    pub(super) fn get_full(&self, value: &str) -> Option<(usize, &str)> {
        self.get_index_of(value)
            .and_then(|i| self.get_index(i).map(|entry| (i, entry)))
    }

    #[doc(hidden)]
    pub fn insert(&mut self, value: String) -> Option<String> {
        self.insert_full(value).1
    }

    fn insert_full(&mut self, value: String) -> (usize, Option<String>) {
        match self.get_index_of(&value) {
            Some(i) => {
                let entry = self.entries[i].replace(value);
                (i, entry)
            }
            None => {
                let i = self.push(value);
                (i, None)
            }
        }
    }

    pub(super) fn insert_at(&mut self, i: usize, value: String) -> Option<String> {
        if i >= self.entries.len() {
            self.entries.resize(i + 1, None);
        }

        self.indices.insert(value.clone(), i);
        self.entries[i].replace(value)
    }

    fn push(&mut self, value: String) -> usize {
        let i = self.entries.len();

        self.indices.insert(value.clone(), i);
        self.entries.push(Some(value));

        i
    }
}
