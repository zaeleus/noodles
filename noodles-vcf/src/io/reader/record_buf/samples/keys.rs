use std::{error, fmt};

use crate::{Header, io::reader::record_buf::MISSING, variant::record_buf::samples::Keys};

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

pub(super) fn parse_keys(header: &Header, s: &str, keys: &mut Keys) -> Result<(), ParseError> {
    const DELIMITER: char = ':';

    if s.is_empty() {
        return Err(ParseError::Empty);
    } else if s == MISSING {
        return Ok(());
    }

    for raw_key in s.split(DELIMITER) {
        let key = match header.formats().get_full(raw_key) {
            Some((_, k, _)) => k.clone(),
            None => raw_key.into(),
        };

        if let Some(key) = keys.as_mut().replace(key) {
            return Err(ParseError::DuplicateKey(key));
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

        let header = Header::default();
        let mut keys = Keys::default();

        keys.as_mut().clear();
        parse_keys(&header, ".", &mut keys)?;
        assert_eq!(keys, Keys::default());

        keys.as_mut().clear();
        parse_keys(&header, "GT", &mut keys)?;
        let expected = [String::from(key::GENOTYPE)].into_iter().collect();
        assert_eq!(keys, expected);

        keys.as_mut().clear();
        parse_keys(&header, "GQ", &mut keys)?;
        let expected = [String::from(key::CONDITIONAL_GENOTYPE_QUALITY)]
            .into_iter()
            .collect();
        assert_eq!(keys, expected);

        keys.as_mut().clear();
        parse_keys(&header, "GT:GQ", &mut keys)?;
        let expected = [
            String::from(key::GENOTYPE),
            String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
        ]
        .into_iter()
        .collect();
        assert_eq!(keys, expected);

        keys.as_mut().clear();
        assert_eq!(parse_keys(&header, "", &mut keys), Err(ParseError::Empty));

        keys.as_mut().clear();
        assert_eq!(
            parse_keys(&header, "GT:GT", &mut keys),
            Err(ParseError::DuplicateKey(String::from(key::GENOTYPE)))
        );

        Ok(())
    }
}
