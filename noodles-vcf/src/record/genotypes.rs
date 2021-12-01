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
pub struct Genotypes {
    keys: Keys,
    genotypes: Vec<Genotype>,
}

impl Genotypes {
    /// Creates VCF record genotypes.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::{genotypes::Keys, Genotypes};
    /// let genotypes = Genotypes::new(Keys::default(), Vec::new());
    /// ```
    pub fn new(keys: Keys, genotypes: Vec<Genotype>) -> Self {
        Self { keys, genotypes }
    }

    /// Returns the genotypes keys.
    pub fn keys(&self) -> &Keys {
        &self.keys
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
        &self.genotypes
    }
}

impl DerefMut for Genotypes {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.genotypes
    }
}

impl fmt::Display for Genotypes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\t", self.keys())?;

        for (i, genotype) in self.genotypes.iter().enumerate() {
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
        let keys = "GT:GQ".parse()?;
        let genotypes = vec![
            Genotype::from_str_format("0|0:7", &keys)?,
            Genotype::from_str_format("./.:20", &keys)?,
            Genotype::from_str_format("1/1:1", &keys)?,
            Genotype::from_str_format(".", &keys)?,
        ];
        let genotypes = Genotypes::new(keys, genotypes);

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

        let genotypes = Genotypes::new(
            Keys::try_from(vec![Key::Genotype, Key::ConditionalGenotypeQuality])?,
            vec![Genotype::try_from(vec![
                Field::new(Key::Genotype, Some(Value::String(String::from("0|0")))),
                Field::new(Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
            ])?],
        );

        assert_eq!(genotypes.to_string(), "GT:GQ\t0|0:13");

        Ok(())
    }
}
