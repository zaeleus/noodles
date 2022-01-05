//! VCF record genotypes and fields.

pub mod genotype;
pub mod keys;

pub use self::{genotype::Genotype, keys::Keys};

use std::{
    error, fmt,
    ops::{Deref, DerefMut},
    str::FromStr,
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

    /// Returns a mutable reference to the genotypes keys.
    pub fn keys_mut(&mut self) -> &mut Keys {
        &mut self.keys
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

/// An error returned when raw VCF record genotypes fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
    /// A key is invalid.
    InvalidKeys(keys::ParseError),
    /// A genotype is invalid.
    InvalidGenotype(genotype::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidKeys(e) => write!(f, "invalid keys: {}", e),
            Self::InvalidGenotype(e) => write!(f, "invalid genotype: {}", e),
        }
    }
}

impl FromStr for Genotypes {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use super::FIELD_DELIMITER;

        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let (format, t) = s.split_once(FIELD_DELIMITER).ok_or(ParseError::Invalid)?;

        let keys = format.parse().map_err(ParseError::InvalidKeys)?;

        let genotypes = t
            .split(FIELD_DELIMITER)
            .map(|t| Genotype::from_str_format(t, &keys))
            .collect::<Result<_, _>>()
            .map_err(ParseError::InvalidGenotype)?;

        Ok(Self::new(keys, genotypes))
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

    #[test]
    fn test_from_str() -> Result<(), Box<dyn std::error::Error>> {
        use genotype::{
            field::{Key, Value},
            Field,
        };

        let expected = Genotypes::new(
            Keys::try_from(vec![Key::Genotype, Key::ConditionalGenotypeQuality])?,
            vec![Genotype::try_from(vec![
                Field::new(Key::Genotype, Some(Value::String(String::from("0|0")))),
                Field::new(Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
            ])?],
        );
        assert_eq!("GT:GQ\t0|0:13".parse(), Ok(expected));

        assert_eq!("".parse::<Genotypes>(), Err(ParseError::Empty));
        assert_eq!("GT:GQ".parse::<Genotypes>(), Err(ParseError::Invalid));
        assert!(matches!(
            "\t".parse::<Genotypes>(),
            Err(ParseError::InvalidKeys(_))
        ));
        assert!(matches!(
            "GQ\tndls".parse::<Genotypes>(),
            Err(ParseError::InvalidGenotype(_))
        ));

        Ok(())
    }
}
