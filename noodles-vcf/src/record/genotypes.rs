//! VCF record genotypes and fields.

pub mod genotype;
pub mod keys;

pub use self::{genotype::Genotype, keys::Keys};

use std::{
    fmt,
    ops::{Deref, DerefMut},
};

use self::genotype::field;

/// VCF record genotypes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Genotypes(Vec<Genotype>);

impl Genotypes {
    /// Creates VCF record genotypes.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::Genotypes;
    /// let genotypes = Genotypes::new(Vec::new());
    /// ```
    pub fn new(genotypes: Vec<Genotype>) -> Self {
        Self(genotypes)
    }

    /// Returns the VCF record genotype value.
    pub fn genotypes(
        &self,
    ) -> Result<Vec<Option<field::value::Genotype>>, genotype::GenotypeError> {
        self.iter().map(|g| g.genotype().transpose()).collect()
    }
}

impl Deref for Genotypes {
    type Target = Vec<Genotype>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Genotypes {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl fmt::Display for Genotypes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, genotype) in self.0.iter().enumerate() {
            if i > 0 {
                write!(f, "\t")?;
            }

            write!(f, "{}", genotype)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genotypes() -> Result<(), Box<dyn std::error::Error>> {
        let format = "GT:GQ".parse()?;

        let genotypes = Genotypes::new(vec![
            Genotype::from_str_format("0|0:7", &format)?,
            Genotype::from_str_format("./.:20", &format)?,
            Genotype::from_str_format("1/1:1", &format)?,
            Genotype::from_str_format(".", &format)?,
        ]);

        let actual = genotypes.genotypes();
        let expected = Ok(vec![
            Some("0|0".parse()?),
            Some("./.".parse()?),
            Some("1/1".parse()?),
            None,
        ]);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_fmt() -> Result<(), Box<dyn std::error::Error>> {
        use genotype::{
            field::{Key, Value},
            Field,
        };

        let genotypes = Genotypes::new(vec![Genotype::try_from(vec![
            Field::new(Key::Genotype, Some(Value::String(String::from("0|0")))),
            Field::new(Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
        ])?]);

        assert_eq!(genotypes.to_string(), "0|0:13");

        Ok(())
    }
}
