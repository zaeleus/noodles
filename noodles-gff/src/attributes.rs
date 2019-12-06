use std::{iter::IntoIterator, str::Split};

const DELIMITER: char = ';';

pub struct Attributes<'a> {
    attrs: &'a str,
}

impl<'a> Attributes<'a> {
    /// Wraps a GFFv2 attribute string for parsing.
    ///
    /// # Example
    ///
    /// ```
    /// use noodles_gff::Attributes;
    ///
    /// let data = r#"gene_name "DDX11L1"; level 2;"#;
    /// let _attributes = Attributes::new(data);
    /// ```
    pub fn new(attrs: &str) -> Attributes {
        Attributes { attrs }
    }

    /// Returns the value corresponding to the key.
    ///
    /// This uses a linear search, which is only really useful for a single
    /// lookup. For multiple lookups, consider converting the attribute pairs
    /// to a map first.
    ///
    /// # Example
    ///
    /// ```
    /// use noodles_gff::Attributes;
    ///
    /// let data = r#"gene_name "DDX11L1"; level 2;"#;
    /// let attributes = Attributes::new(data);
    ///
    /// assert_eq!(attributes.get("gene_name"), Some("DDX11L1"));
    /// assert_eq!(attributes.get("level"), Some("2"));
    /// assert_eq!(attributes.get("gene_type"), None);
    /// ```
    pub fn get(&self, key: &str) -> Option<&str> {
        self.iter().find(|&(k, _)| k == key).map(|(_, v)| v)
    }

    /// Returns an iterator that parses over all attribute pairs.
    ///
    /// # Example
    ///
    /// ```
    /// use noodles_gff::Attributes;
    ///
    /// let data = r#"gene_name "DDX11L1"; level 2;"#;
    /// let attributes = Attributes::new(data);
    /// let mut it = attributes.iter();
    ///
    /// assert_eq!(it.next(), Some(("gene_name", "DDX11L1")));
    /// assert_eq!(it.next(), Some(("level", "2")));
    /// assert_eq!(it.next(), None)
    /// ```
    pub fn iter(&self) -> AttributesIter {
        self.into_iter()
    }
}

impl<'a> IntoIterator for &'a Attributes<'a> {
    type Item = (&'a str, &'a str);
    type IntoIter = AttributesIter<'a>;

    fn into_iter(self) -> AttributesIter<'a> {
        AttributesIter {
            split: self.attrs.split(DELIMITER),
        }
    }
}

pub struct AttributesIter<'a> {
    split: Split<'a, char>,
}

impl<'a> Iterator for AttributesIter<'a> {
    type Item = (&'a str, &'a str);

    fn next(&mut self) -> Option<Self::Item> {
        self.split.next().map(str::trim_start).and_then(|p| {
            let mut pieces = p.splitn(2, ' ');

            if let Some(key) = pieces.next() {
                if let Some(value) = pieces.next() {
                    return Some((key, trim_quotes(value)));
                }
            }

            None
        })
    }
}

fn trim_quotes(s: &str) -> &str {
    s.trim_matches('"')
}

#[cfg(test)]
mod attributes_iter_tests {
    use super::Attributes;

    #[test]
    fn test_next() {
        let data = r#"gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";"#;
        let attributes = Attributes::new(data);
        let mut it = attributes.iter();

        assert_eq!(it.next(), Some(("gene_id", "ENSG00000223972.5")));
        assert_eq!(
            it.next(),
            Some(("gene_type", "transcribed_unprocessed_pseudogene"))
        );
        assert_eq!(it.next(), Some(("gene_name", "DDX11L1")));
        assert_eq!(it.next(), Some(("level", "2")));
        assert_eq!(it.next(), Some(("havana_gene", "OTTHUMG00000000961.2")));
        assert_eq!(it.next(), None);
    }
}

#[cfg(test)]
mod tests {
    use super::trim_quotes;

    #[test]
    fn test_trim_quotes() {
        assert_eq!(trim_quotes("DDX11L1"), "DDX11L1");
        assert_eq!(trim_quotes(r#""DDX11L1""#), "DDX11L1");
        assert_eq!(trim_quotes(""), "");
    }
}
