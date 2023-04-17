//! VCF record genotypes and fields.

pub mod keys;
pub mod sample;

pub use self::{keys::Keys, sample::Sample};

use std::{
    error,
    fmt::{self, Write},
    str::FromStr,
};

use self::sample::Value;
use super::FIELD_DELIMITER;
use crate::{
    header::{
        record::value::{map::Format, Map},
        Formats,
    },
    Header,
};

/// VCF record genotypes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Genotypes {
    pub(crate) keys: Keys,
    pub(crate) values: Vec<Vec<Option<Value>>>,
}

impl Genotypes {
    /// Parses VCF record genotypes with a VCF header as context.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::record::value::{map::Format, Map},
    ///     record::{
    ///         genotypes::{keys::key, sample::Value, Keys},
    ///         Genotypes,
    ///     },
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_format(key::GENOTYPE, Map::<Format>::from(&key::GENOTYPE))
    ///     .build();
    ///
    /// let actual = Genotypes::parse("GT:GQ\t0|0:13", &header)?;
    ///
    /// let expected = Genotypes::new(
    ///     Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?,
    ///     vec![vec![
    ///         Some(Value::String(String::from("0|0"))),
    ///         Some(Value::Integer(13)),
    ///     ]],
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
    /// let keys = Keys::try_from(vec![key::GENOTYPE])?;
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
    /// let keys = Keys::try_from(vec![key::GENOTYPE])?;
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

impl fmt::Display for Genotypes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.keys(), FIELD_DELIMITER)?;

        for (i, sample) in self.values().enumerate() {
            if i > 0 {
                f.write_char(FIELD_DELIMITER)?;
            }

            for (j, value) in sample.values().iter().enumerate() {
                if j > 0 {
                    ':'.fmt(f)?;
                }

                if let Some(v) = value {
                    write!(f, "{v}")?;
                } else {
                    '.'.fmt(f)?;
                }
            }
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

    let (format, rest) = s.split_once(FIELD_DELIMITER).ok_or(ParseError::Invalid)?;

    let keys = format.parse().map_err(ParseError::InvalidKeys)?;

    let values = rest
        .split(FIELD_DELIMITER)
        .map(|t| parse_values(t, header.formats(), &keys))
        .collect::<Result<_, _>>()
        .map_err(ParseError::InvalidValues)?;

    Ok(Genotypes::new(keys, values))
}

fn parse_values(
    s: &str,
    formats: &Formats,
    keys: &Keys,
) -> Result<Vec<Option<Value>>, sample::ParseError> {
    if s.is_empty() {
        return Err(sample::ParseError::Empty);
    } else if s == "." {
        return Ok(Vec::new());
    }

    let mut values = Vec::with_capacity(keys.len());
    let mut raw_values = s.split(':');

    for (key, raw_value) in keys.iter().zip(&mut raw_values) {
        let value = if let Some(format) = formats.get(key) {
            parse_value(format, raw_value).map_err(sample::ParseError::InvalidValue)?
        } else {
            let format = Map::<Format>::from(key);
            parse_value(&format, raw_value).map_err(sample::ParseError::InvalidValue)?
        };

        values.push(value);
    }

    if raw_values.next().is_some() {
        Err(sample::ParseError::UnexpectedValue)
    } else {
        Ok(values)
    }
}

fn parse_value(format: &Map<Format>, s: &str) -> Result<Option<Value>, sample::value::ParseError> {
    match s {
        "." => Ok(None),
        _ => sample::Value::from_str_format(s, format).map(Some),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::genotypes::keys::key;

    #[test]
    fn test_genotypes() -> Result<(), Box<dyn std::error::Error>> {
        use crate::header::record::value::{map::Format, Map};

        let header = crate::Header::builder()
            .add_format(key::GENOTYPE, Map::<Format>::from(&key::GENOTYPE))
            .add_format(
                key::CONDITIONAL_GENOTYPE_QUALITY,
                Map::<Format>::from(&key::CONDITIONAL_GENOTYPE_QUALITY),
            )
            .build();

        let genotypes = Genotypes::parse("GT:GQ\t0|0:7\t./.:20\t1/1:1\t.", &header)?;

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
        let genotypes = Genotypes::new(
            Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?,
            vec![vec![
                Some(Value::String(String::from("0|0"))),
                Some(Value::Integer(13)),
            ]],
        );

        assert_eq!(genotypes.to_string(), "GT:GQ\t0|0:13");

        Ok(())
    }

    #[test]
    fn test_from_str() -> Result<(), super::keys::TryFromKeyVectorError> {
        let expected = Genotypes::new(
            Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?,
            vec![vec![
                Some(Value::String(String::from("0|0"))),
                Some(Value::Integer(13)),
            ]],
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
            Err(ParseError::InvalidValues(_))
        ));

        Ok(())
    }
}
