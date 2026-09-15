use std::{error, fmt};

use crate::{io::reader::record_buf::MISSING, variant::record_buf::samples::Keys};

/// An error when raw VCF record genotypes keys fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A key is duplicated.
    DuplicateKey(String),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::DuplicateKey(key) => write!(f, "duplicate key: {key}"),
        }
    }
}

pub(super) fn parse_keys(s: &str, keys: &mut Keys) -> Result<(), ParseError> {
    const DELIMITER: char = ':';

    if s.is_empty() {
        return Err(ParseError::Empty);
    } else if s == MISSING {
        return Ok(());
    }

    for key in s.split(DELIMITER) {
        if let Some(k) = keys.as_mut().replace(key.into()) {
            return Err(ParseError::DuplicateKey(k));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_keys() -> Result<(), Box<dyn std::error::Error>> {
        use crate::variant::record::samples::keys::key;

        let mut keys = Keys::default();

        keys.as_mut().clear();
        parse_keys(".", &mut keys)?;
        assert_eq!(keys, Keys::default());

        keys.as_mut().clear();
        parse_keys("GT", &mut keys)?;
        let expected = [String::from(key::GENOTYPE)].into_iter().collect();
        assert_eq!(keys, expected);

        keys.as_mut().clear();
        parse_keys("GQ", &mut keys)?;
        let expected = [String::from(key::CONDITIONAL_GENOTYPE_QUALITY)]
            .into_iter()
            .collect();
        assert_eq!(keys, expected);

        keys.as_mut().clear();
        parse_keys("GT:GQ", &mut keys)?;
        let expected = [
            String::from(key::GENOTYPE),
            String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
        ]
        .into_iter()
        .collect();
        assert_eq!(keys, expected);

        keys.as_mut().clear();
        assert_eq!(parse_keys("", &mut keys), Err(ParseError::Empty));

        keys.as_mut().clear();
        assert_eq!(
            parse_keys("GT:GT", &mut keys),
            Err(ParseError::DuplicateKey(String::from(key::GENOTYPE)))
        );

        Ok(())
    }
}
