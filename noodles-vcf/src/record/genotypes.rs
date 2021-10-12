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
    ) -> Option<Result<Vec<field::value::Genotype>, field::value::genotype::ParseError>> {
        self.iter().map(Genotype::genotype).collect()
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
    use std::convert::TryFrom;

    use super::*;

    #[test]
    fn test_genotypes() -> Result<(), Box<dyn std::error::Error>> {
        let format = "GT:GQ".parse()?;

        let genotypes: Genotypes = vec![
            Genotype::from_str_format("0|0:7", &format)?,
            Genotype::from_str_format("./.:20", &format)?,
            Genotype::from_str_format("1/1:1", &format)?,
        ]
        .into();

        let genotype_values = genotypes.genotypes();

        let expected = vec!["0|0".parse()?, "./.".parse()?, "1/1".parse()?];

        assert_eq!(genotype_values, Some(Ok(expected)));

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
