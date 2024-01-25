mod keys;
mod values;

use std::{error, fmt};

use self::{keys::parse_keys, values::parse_values};
use super::next_field;
use crate::{record::Genotypes, Header};

/// An error when raw VCF record genotypes fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is unexpected.
    ///
    /// The header has no samples, but the record has unexpected data.
    UnexpectedInput,
    /// The keys are invalid.
    InvalidKeys(keys::ParseError),
    /// A list of sample values is invalid.
    InvalidValues(values::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::UnexpectedInput => None,
            Self::InvalidKeys(e) => Some(e),
            Self::InvalidValues(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::UnexpectedInput => write!(f, "unexpected input"),
            ParseError::InvalidKeys(_) => write!(f, "invalid keys"),
            ParseError::InvalidValues(_) => write!(f, "invalid values"),
        }
    }
}

pub(super) fn parse_genotypes(
    header: &Header,
    mut s: &str,
    genotypes: &mut Genotypes,
) -> Result<(), ParseError> {
    genotypes.keys.clear();

    let sample_count = header.sample_names().len();

    if sample_count == 0 {
        genotypes.values.clear();

        if s.is_empty() {
            return Ok(());
        } else {
            return Err(ParseError::UnexpectedInput);
        }
    }

    for values in &mut genotypes.values {
        values.clear();
    }

    let field = next_field(&mut s);
    parse_keys(header, field, &mut genotypes.keys).map_err(ParseError::InvalidKeys)?;

    genotypes.values.resize(sample_count, Vec::new());

    for values in &mut genotypes.values {
        let field = next_field(&mut s);
        parse_values(header, &genotypes.keys, field, values).map_err(ParseError::InvalidValues)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_genotypes() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::genotypes::{keys::key, sample::Value, Keys};

        let mut genotypes = Genotypes::default();

        let header = Header::default();
        parse_genotypes(&header, "", &mut genotypes)?;
        assert!(genotypes.is_empty());

        let header = Header::builder().add_sample_name("sample0").build();
        parse_genotypes(&header, "GT\t0|0", &mut genotypes)?;
        let expected = Genotypes::new(
            Keys::try_from(vec![key::GENOTYPE])?,
            vec![vec![Some(Value::from("0|0"))]],
        );
        assert_eq!(genotypes, expected);

        let header = Header::builder()
            .add_sample_name("sample0")
            .add_sample_name("sample1")
            .build();
        parse_genotypes(&header, "GQ\t8\t13", &mut genotypes)?;
        let expected = Genotypes::new(
            Keys::try_from(vec![key::CONDITIONAL_GENOTYPE_QUALITY])?,
            vec![vec![Some(Value::from(8))], vec![Some(Value::from(13))]],
        );
        assert_eq!(genotypes, expected);

        let header = Header::default();
        assert_eq!(
            parse_genotypes(&header, "GT\t0|0", &mut genotypes),
            Err(ParseError::UnexpectedInput)
        );

        let header = Header::builder().add_sample_name("sample0").build();

        assert!(matches!(
            parse_genotypes(&header, "\t0|0", &mut genotypes),
            Err(ParseError::InvalidKeys(_))
        ));

        assert!(matches!(
            parse_genotypes(&header, "GT:GQ", &mut genotypes),
            Err(ParseError::InvalidValues(_))
        ));

        assert!(matches!(
            parse_genotypes(&header, "GQ\tndls", &mut genotypes),
            Err(ParseError::InvalidValues(_))
        ));

        Ok(())
    }
}
