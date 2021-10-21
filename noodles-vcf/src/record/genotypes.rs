//! VCF record genotypes and fields.

pub mod genotype;

pub use self::genotype::Genotype;

use std::{
    fmt,
    ops::{Deref, DerefMut},
};

use self::genotype::field;

/// VCF record genotypes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Genotypes(Vec<Genotype>);

impl Genotypes {
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
        for (i, field) in self.0.iter().enumerate() {
            if i > 0 {
                write!(f, "\t")?;
            }

            write!(f, "{}", field)?;
        }

        Ok(())
    }
}

impl From<Vec<Genotype>> for Genotypes {
    fn from(genotypes: Vec<Genotype>) -> Self {
        Self(genotypes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genotypes() -> Result<(), Box<dyn std::error::Error>> {
        let format = "GT:GQ".parse()?;

        let genotypes = Genotypes::from(vec![
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

        let genotypes = Genotypes::from(vec![Genotype::try_from(vec![
            Field::new(Key::Genotype, Some(Value::String(String::from("0|0")))),
            Field::new(Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
        ])?]);

        assert_eq!(genotypes.to_string(), "0|0:13");

        Ok(())
    }
}
