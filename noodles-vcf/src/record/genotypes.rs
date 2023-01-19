//! VCF record genotypes and fields.

pub mod genotype;
pub mod keys;

pub use self::{genotype::Genotype, keys::Keys};

use std::{
    error,
    fmt::{self, Write},
    ops::{Deref, DerefMut},
    str::FromStr,
};

use self::genotype::field;
use super::FIELD_DELIMITER;
use crate::Header;

/// VCF record genotypes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Genotypes {
    keys: Keys,
    genotypes: Vec<Genotype>,
}

impl Genotypes {
    /// Parses VCF record genotypes with a VCF header as context.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::{format::Key, record::value::{map::Format, Map}},
    ///     record::{
    ///         genotypes::{genotype::field::Value, Keys},
    ///         Genotypes,
    ///     },
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_format(Key::Genotype, Map::<Format>::from(&Key::Genotype))
    ///     .build();
    ///
    /// let actual = Genotypes::parse("GT:GQ\t0|0:13", &header)?;
    ///
    /// let expected = Genotypes::new(
    ///     Keys::try_from(vec![Key::Genotype, Key::ConditionalGenotypeQuality])?,
    ///     vec![[
    ///         (Key::Genotype, Some(Value::String(String::from("0|0")))),
    ///         (Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
    ///     ].into_iter().collect()],
    /// );
    ///
    /// assert_eq!(actual, expected);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn parse(s: &str, header: &Header) -> Result<Genotypes, ParseError> {
        parse(s, header)
    }

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
        self.genotypes.is_empty()
    }

    /// Returns the genotypes keys.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{header::format::Key, record::{genotypes::Keys, Genotypes}};
    ///
    /// let genotypes = Genotypes::default();
    /// assert!(genotypes.keys().is_empty());
    ///
    /// let keys = Keys::try_from(vec![Key::Genotype])?;
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
    /// use noodles_vcf::{header::format::Key, record::{genotypes::Keys, Genotypes}};
    ///
    /// let keys = Keys::try_from(vec![Key::Genotype])?;
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
        write!(f, "{}{}", self.keys(), FIELD_DELIMITER)?;

        for (i, genotype) in self.genotypes.iter().enumerate() {
            if i > 0 {
                f.write_char(FIELD_DELIMITER)?;
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

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKeys(e) => Some(e),
            Self::InvalidGenotype(e) => Some(e),
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
            Self::InvalidGenotype(_) => f.write_str("invalid genotype"),
        }
    }
}

impl FromStr for Genotypes {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse(s, &Header::default())
    }
}

fn parse(s: &str, header: &Header) -> Result<Genotypes, ParseError> {
    if s.is_empty() {
        return Err(ParseError::Empty);
    }

    let (format, t) = s.split_once(FIELD_DELIMITER).ok_or(ParseError::Invalid)?;

    let keys = format.parse().map_err(ParseError::InvalidKeys)?;

    let genotypes = t
        .split(FIELD_DELIMITER)
        .map(|t| Genotype::parse(t, header.formats(), &keys))
        .collect::<Result<_, _>>()
        .map_err(ParseError::InvalidGenotype)?;

    Ok(Genotypes::new(keys, genotypes))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genotypes() -> Result<(), Box<dyn std::error::Error>> {
        use crate::header::{
            format::Key,
            record::value::{map::Format, Map},
        };

        let header = crate::Header::builder()
            .add_format(Key::Genotype, Map::<Format>::from(&Key::Genotype))
            .add_format(
                Key::ConditionalGenotypeQuality,
                Map::<Format>::from(&Key::ConditionalGenotypeQuality),
            )
            .build();

        let keys = "GT:GQ".parse()?;
        let genotypes = vec![
            Genotype::parse("0|0:7", header.formats(), &keys)?,
            Genotype::parse("./.:20", header.formats(), &keys)?,
            Genotype::parse("1/1:1", header.formats(), &keys)?,
            Genotype::parse(".", header.formats(), &keys)?,
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
    fn test_fmt() -> Result<(), super::keys::TryFromKeyVectorError> {
        use self::genotype::field::Value;
        use crate::header::format::Key;

        let genotypes = Genotypes::new(
            Keys::try_from(vec![Key::Genotype, Key::ConditionalGenotypeQuality])?,
            vec![[
                (Key::Genotype, Some(Value::String(String::from("0|0")))),
                (Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
            ]
            .into_iter()
            .collect()],
        );

        assert_eq!(genotypes.to_string(), "GT:GQ\t0|0:13");

        Ok(())
    }

    #[test]
    fn test_from_str() -> Result<(), super::keys::TryFromKeyVectorError> {
        use super::genotype::field::Value;
        use crate::header::format::Key;

        let expected = Genotypes::new(
            Keys::try_from(vec![Key::Genotype, Key::ConditionalGenotypeQuality])?,
            vec![[
                (Key::Genotype, Some(Value::String(String::from("0|0")))),
                (Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
            ]
            .into_iter()
            .collect()],
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
