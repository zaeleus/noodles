//! VCF header meta record.

use std::fmt;

use super::record;

/// A VCF header meta record (`META`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Meta {
    id: String,
    values: Vec<String>,
}

impl Meta {
    /// Creates a VCF header meta record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Meta;
    ///
    /// let meta = Meta::new(
    ///     String::from("Assay"),
    ///     vec![String::from("WholeGenome"), String::from("Exome")],
    /// );
    /// ```
    pub fn new(id: String, values: Vec<String>) -> Self {
        Self { id, values }
    }

    /// Returns the ID of the meta record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Meta;
    ///
    /// let meta = Meta::new(
    ///     String::from("Assay"),
    ///     vec![String::from("WholeGenome"), String::from("Exome")],
    /// );
    ///
    /// assert_eq!(meta.id(), "Assay");
    /// ```
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns the values of the meta record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Meta;
    ///
    /// let meta = Meta::new(
    ///     String::from("Assay"),
    ///     vec![String::from("WholeGenome"), String::from("Exome")],
    /// );
    ///
    /// assert_eq!(meta.values(), ["WholeGenome", "Exome"]);
    /// ```
    pub fn values(&self) -> &[String] {
        &self.values
    }
}

impl fmt::Display for Meta {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::Key::Meta.as_ref())?;
        f.write_str("=<")?;

        write!(f, "ID={}", self.id)?;
        f.write_str(",Type=String")?;
        f.write_str(",Number=.")?;

        f.write_str(",Values=[")?;

        for (i, value) in self.values().iter().enumerate() {
            if i > 0 {
                f.write_str(", ")?;
            }

            f.write_str(value)?;
        }

        f.write_str("]")?;

        f.write_str(">")?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let meta = Meta::new(
            String::from("Assay"),
            vec![String::from("WholeGenome"), String::from("Exome")],
        );

        let expected = r#"##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>"#;
        assert_eq!(meta.to_string(), expected);
    }
}
