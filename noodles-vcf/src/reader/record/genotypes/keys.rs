use std::{error, fmt};

use crate::{
    reader::record::MISSING,
    record::genotypes::{
        keys::{key, Key},
        Keys,
    },
    Header,
};

/// An error when raw VCF record genotypes keys fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A key is invalid.
    InvalidKey(key::ParseError),
    /// The genotype key (`GT`) position is invalid.
    ///
    /// The genotype key must be first, if present.
    InvalidGenotypeKeyPosition,
    /// A key is duplicated.
    DuplicateKey(Key),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKey(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::InvalidKey(_) => write!(f, "invalid key"),
            Self::InvalidGenotypeKeyPosition => write!(f, "invalid genotype key position"),
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

    let mut gt_position = None;

    for (i, raw_key) in s.split(DELIMITER).enumerate() {
        let key = match header.formats().get_full(raw_key) {
            Some((_, k, _)) => k.clone(),
            None => {
                Key::try_from((header.file_format(), raw_key)).map_err(ParseError::InvalidKey)?
            }
        };

        if key == key::GENOTYPE {
            gt_position = Some(i);
        }

        if let Some(key) = keys.replace(key) {
            return Err(ParseError::DuplicateKey(key));
        }
    }

    if let Some(i) = gt_position {
        if i != 0 {
            return Err(ParseError::InvalidGenotypeKeyPosition);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_keys() -> Result<(), Box<dyn std::error::Error>> {
        let header = Header::default();
        let mut keys = Keys::default();

        keys.clear();
        parse_keys(&header, ".", &mut keys)?;
        assert_eq!(keys, Keys::default());

        keys.clear();
        parse_keys(&header, "GT", &mut keys)?;
        let expected = Keys::try_from(vec![key::GENOTYPE])?;
        assert_eq!(keys, expected);

        keys.clear();
        parse_keys(&header, "GQ", &mut keys)?;
        let expected = Keys::try_from(vec![key::CONDITIONAL_GENOTYPE_QUALITY])?;
        assert_eq!(keys, expected);

        keys.clear();
        parse_keys(&header, "GT:GQ", &mut keys)?;
        let expected = Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?;
        assert_eq!(keys, expected);

        keys.clear();
        assert_eq!(parse_keys(&header, "", &mut keys), Err(ParseError::Empty));

        keys.clear();
        assert_eq!(
            parse_keys(&header, "GQ:GT", &mut keys),
            Err(ParseError::InvalidGenotypeKeyPosition)
        );

        keys.clear();
        assert_eq!(
            parse_keys(&header, "GT:GT", &mut keys),
            Err(ParseError::DuplicateKey(key::GENOTYPE))
        );

        Ok(())
    }
}
