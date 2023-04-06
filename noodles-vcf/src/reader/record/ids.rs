mod id;

use std::{error, fmt};

use noodles_core as core;

use self::id::parse_id;
use crate::record::Ids;

/// An error when raw VCF record IDs fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// An ID is invalid.
    InvalidId(id::ParseError),
    /// An ID is duplicated.
    DuplicateId,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::InvalidId(_) => write!(f, "invalid ID"),
            Self::DuplicateId => write!(f, "duplicate ID"),
        }
    }
}

impl From<ParseError> for core::Error {
    fn from(e: ParseError) -> Self {
        Self::new(core::error::Kind::Parse, e)
    }
}

pub(super) fn parse_ids(s: &str, ids: &mut Ids) -> Result<(), ParseError> {
    const DELIMITER: char = ';';

    if s.is_empty() {
        return Err(ParseError::Empty);
    }

    for raw_id in s.split(DELIMITER) {
        let id = parse_id(raw_id).map_err(ParseError::InvalidId)?;

        if !ids.insert(id) {
            return Err(ParseError::DuplicateId);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_ids() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::ids::Id;

        let id0: Id = "nd0".parse()?;
        let id1: Id = "nd1".parse()?;

        let mut ids = Ids::default();

        ids.clear();
        parse_ids("nd0", &mut ids)?;
        let expected = [id0.clone()].into_iter().collect();
        assert_eq!(ids, expected);

        ids.clear();
        parse_ids("nd0;nd1", &mut ids)?;
        let expected = [id0, id1].into_iter().collect();
        assert_eq!(ids, expected);

        ids.clear();
        assert_eq!(parse_ids("", &mut ids), Err(ParseError::Empty));

        ids.clear();
        assert!(matches!(
            parse_ids("nd 0", &mut ids),
            Err(ParseError::InvalidId(_))
        ));

        ids.clear();
        assert!(matches!(
            parse_ids(";nd0", &mut ids),
            Err(ParseError::InvalidId(_))
        ));

        ids.clear();
        assert!(matches!(
            parse_ids("nd0;;nd1", &mut ids),
            Err(ParseError::InvalidId(_))
        ));

        ids.clear();
        assert!(matches!(
            parse_ids("nd0;", &mut ids),
            Err(ParseError::InvalidId(_))
        ));

        ids.clear();
        assert_eq!(parse_ids("nd0;nd0", &mut ids), Err(ParseError::DuplicateId));

        Ok(())
    }
}
