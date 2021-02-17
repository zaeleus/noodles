//! VCF header sample record.

use std::{collections::HashMap, fmt};

use super::record;

/// A VCF header sample record (`SAMPLE`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Sample {
    id: String,
    fields: HashMap<String, String>,
}

impl Sample {
    /// Creates a VCF header sample record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Sample;
    /// let sample = Sample::new(String::from("sample0"), Default::default());
    /// ```
    pub fn new(id: String, fields: HashMap<String, String>) -> Self {
        Self { id, fields }
    }

    /// Returns the ID of the sample record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Sample;
    /// let sample = Sample::new(String::from("sample0"), Default::default());
    /// assert_eq!(sample.id(), "sample0");
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
    /// use noodles_vcf::header::Sample;
    /// let sample = Sample::new(String::from("sample0"), Default::default());
    /// assert!(sample.fields().is_empty());
    /// ```
    pub fn fields(&self) -> &HashMap<String, String> {
        &self.fields
    }
}

impl fmt::Display for Sample {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::Key::Sample.as_ref())?;
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
        let sample = Sample::new(String::from("sample0"), HashMap::new());
        assert_eq!(sample.to_string(), "##SAMPLE=<ID=sample0>");

        let mut fields = HashMap::new();
        fields.insert(String::from("Assay"), String::from("WholeGenome"));
        let sample = Sample::new(String::from("sample0"), fields);
        assert_eq!(
            sample.to_string(),
            "##SAMPLE=<ID=sample0,Assay=WholeGenome>"
        );
    }
}
