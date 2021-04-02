//! VCF header pedigree record.

use std::fmt;

use indexmap::IndexMap;

use super::record;

/// A VCF header pedigree record (`PEDIGREE`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Pedigree {
    id: String,
    fields: IndexMap<String, String>,
}

impl Pedigree {
    /// Creates a VCF header pedigree record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Pedigree;
    /// let pedigree = Pedigree::new(String::from("cid"), Default::default());
    /// ```
    pub fn new(id: String, fields: IndexMap<String, String>) -> Self {
        Self { id, fields }
    }

    /// Returns the ID of the sample record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Pedigree;
    /// let pedigree = Pedigree::new(String::from("cid"), Default::default());
    /// assert_eq!(pedigree.id(), "cid");
    /// ```
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns the extra fields in the record.
    ///
    /// This includes fields other than `ID`.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Pedigree;
    /// let pedigree = Pedigree::new(String::from("cid"), Default::default());
    /// assert!(pedigree.fields().is_empty());
    /// ```
    pub fn fields(&self) -> &IndexMap<String, String> {
        &self.fields
    }
}

impl fmt::Display for Pedigree {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::Key::Pedigree.as_ref())?;
        f.write_str("=<")?;

        write!(f, "ID={}", self.id())?;

        for (key, value) in &self.fields {
            write!(f, r#",{}={}"#, key, value)?;
        }

        f.write_str(">")?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let sample = Pedigree::new(String::from("cid"), IndexMap::new());
        assert_eq!(sample.to_string(), "##PEDIGREE=<ID=cid>");

        let mut fields = IndexMap::new();
        fields.insert(String::from("Father"), String::from("fid"));
        fields.insert(String::from("Mother"), String::from("mid"));
        let sample = Pedigree::new(String::from("cid"), fields);
        assert_eq!(
            sample.to_string(),
            "##PEDIGREE=<ID=cid,Father=fid,Mother=mid>"
        );
    }
}
