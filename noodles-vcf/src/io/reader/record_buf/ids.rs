mod id;

use std::{error, fmt};

use self::id::parse_id;
use crate::variant::record_buf::Ids;

/// An error when raw VCF record IDs fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// An ID is invalid.
    InvalidId(id::ParseError),
    /// An ID is duplicated.
    DuplicateId(String),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidId(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::InvalidId(_) => write!(f, "invalid ID"),
            Self::DuplicateId(id) => write!(f, "duplicate ID: {id}"),
        }
    }
}

pub(super) fn parse_ids(s: &str, ids: &mut Ids) -> Result<(), ParseError> {
    const DELIMITER: char = ';';

    if s.is_empty() {
        return Err(ParseError::Empty);
    }

    for raw_id in s.split(DELIMITER) {
        let id = parse_id(raw_id).map_err(ParseError::InvalidId)?;

        if let Some(id) = ids.as_mut().replace(id.into()) {
            return Err(ParseError::DuplicateId(id));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_ids() -> Result<(), Box<dyn std::error::Error>> {
        let id0 = String::from("nd0");
        let id1 = String::from("nd1");

        let mut ids = Ids::default();

        ids.as_mut().clear();
        parse_ids("nd0", &mut ids)?;
        let expected = [id0.clone()].into_iter().collect();
        assert_eq!(ids, expected);

        ids.as_mut().clear();
        parse_ids("nd0;nd1", &mut ids)?;
        let expected = [id0.clone(), id1].into_iter().collect();
        assert_eq!(ids, expected);

        ids.as_mut().clear();
        assert_eq!(parse_ids("", &mut ids), Err(ParseError::Empty));

        ids.as_mut().clear();
        assert!(matches!(
            parse_ids(";nd0", &mut ids),
            Err(ParseError::InvalidId(_))
        ));

        ids.as_mut().clear();
        assert!(matches!(
            parse_ids("nd0;;nd1", &mut ids),
            Err(ParseError::InvalidId(_))
        ));

        ids.as_mut().clear();
        assert!(matches!(
            parse_ids("nd0;", &mut ids),
            Err(ParseError::InvalidId(_))
        ));

        ids.as_mut().clear();
        assert_eq!(
            parse_ids("nd0;nd0", &mut ids),
            Err(ParseError::DuplicateId(id0))
        );

        Ok(())
    }
}
