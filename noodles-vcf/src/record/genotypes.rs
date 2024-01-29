//! VCF record genotypes and fields.

pub mod keys;
pub mod sample;

pub use self::{keys::Keys, sample::Sample};

use std::{error, fmt};

use self::sample::Value;

/// VCF record genotypes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Genotypes {
    pub(crate) keys: Keys,
    pub(crate) values: Vec<Vec<Option<Value>>>,
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
    pub fn new(keys: Keys, values: Vec<Vec<Option<Value>>>) -> Self {
        Self { keys, values }
    }

    /// Returns whether there are any samples.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::Genotypes;
    /// let genotypes = Genotypes::default();
    /// assert!(genotypes.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Returns the genotypes keys.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::{genotypes::{keys::key, Keys}, Genotypes};
    ///
    /// let genotypes = Genotypes::default();
    /// assert!(genotypes.keys().is_empty());
    ///
    /// let keys = Keys::try_from(vec![String::from(key::GENOTYPE)])?;
    /// let genotypes = Genotypes::new(keys.clone(), Vec::new());
    /// assert_eq!(genotypes.keys(), &keys);
    /// # Ok::<_, noodles_vcf::record::genotypes::keys::TryFromKeyVectorError>(())
    /// ```
    pub fn keys(&self) -> &Keys {
        &self.keys
    }

    /// Returns a mutable reference to the genotypes keys.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::{genotypes::{keys::key, Keys}, Genotypes};
    ///
    /// let keys = Keys::try_from(vec![String::from(key::GENOTYPE)])?;
    ///
    /// let mut genotypes = Genotypes::default();
    /// *genotypes.keys_mut() = keys.clone();
    ///
    /// assert_eq!(genotypes.keys(), &keys);
    /// # Ok::<_, noodles_vcf::record::genotypes::keys::TryFromKeyVectorError>(())
    /// ```
    pub fn keys_mut(&mut self) -> &mut Keys {
        &mut self.keys
    }

    /// Returns genotypes samples.
    pub fn values(&self) -> impl Iterator<Item = Sample<'_>> {
        self.values
            .iter()
            .map(|values| Sample::new(&self.keys, values))
    }

    /// Returns the genotype values for the sample at the given index.
    pub fn get_index(&self, i: usize) -> Option<Sample<'_>> {
        self.values
            .get(i)
            .map(|values| Sample::new(&self.keys, values))
    }

    /// Returns the VCF record genotype value.
    pub fn genotypes(&self) -> Result<Vec<Option<sample::value::Genotype>>, sample::GenotypeError> {
        self.values()
            .map(|sample| sample.genotype().transpose())
            .collect()
    }
}

/// An error returned when raw VCF record genotypes fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
    /// A key is invalid.
    InvalidKeys(keys::ParseError),
    /// A value is invalid.
    InvalidValues(sample::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKeys(e) => Some(e),
            Self::InvalidValues(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidKeys(_) => f.write_str("invalid keys"),
            Self::InvalidValues(_) => f.write_str("invalid values"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::genotypes::keys::key;

    #[test]
    fn test_genotypes() -> Result<(), Box<dyn std::error::Error>> {
        let keys = Keys::try_from(vec![
            String::from(key::GENOTYPE),
            String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
        ])?;
        let genotypes = Genotypes::new(
            keys,
            vec![
                vec![Some(Value::from("0|0")), Some(Value::from(7))],
                vec![Some(Value::from("./.")), Some(Value::from(20))],
                vec![Some(Value::from("1/1")), Some(Value::from(1))],
                vec![],
            ],
        );

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
}
