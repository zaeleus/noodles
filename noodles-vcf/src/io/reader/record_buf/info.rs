mod field;

use std::{error, fmt};

use self::field::parse_field;
use crate::{Header, variant::record_buf::Info};

/// An error when raw VCF record info fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A field is invalid.
    InvalidField(field::ParseError),
    /// A key is duplicated.
    DuplicateKey(String),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::Empty => write!(f, "empty input"),
            ParseError::InvalidField(e) => {
                write!(f, "invalid field")?;

                if let Some(key) = e.key() {
                    write!(f, ": {key}")?;
                }

                Ok(())
            }
            ParseError::DuplicateKey(key) => write!(f, "duplicate key: {key}"),
        }
    }
}

pub(super) fn parse_info(header: &Header, s: &str, info: &mut Info) -> Result<(), ParseError> {
    use indexmap::map::Entry;

    const DELIMITER: char = ';';

    if s.is_empty() {
        return Err(ParseError::Empty);
    }

    for raw_field in s.split(DELIMITER) {
        let (key, value) = parse_field(header, raw_field).map_err(ParseError::InvalidField)?;

        match info.as_mut().entry(key) {
            Entry::Vacant(entry) => {
                entry.insert(value);
            }
            Entry::Occupied(entry) => {
                let (k, _) = entry.swap_remove_entry();
                return Err(ParseError::DuplicateKey(k));
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_info() -> Result<(), ParseError> {
        use crate::variant::{record::info::field::key, record_buf::info::field::Value};

        let header = Header::default();
        let mut info = Info::default();

        info.clear();
        parse_info(&header, "NS=2", &mut info)?;
        let expected = [(
            String::from(key::SAMPLES_WITH_DATA_COUNT),
            Some(Value::from(2)),
        )]
        .into_iter()
        .collect();
        assert_eq!(info, expected);

        info.clear();
        parse_info(&header, "NS=2;AA=T", &mut info)?;
        let expected = [
            (
                String::from(key::SAMPLES_WITH_DATA_COUNT),
                Some(Value::from(2)),
            ),
            (String::from(key::ANCESTRAL_ALLELE), Some(Value::from("T"))),
        ]
        .into_iter()
        .collect();
        assert_eq!(info, expected);

        assert_eq!(parse_info(&header, "", &mut info), Err(ParseError::Empty));

        assert_eq!(
            parse_info(&header, "NS=2;NS=2", &mut info),
            Err(ParseError::DuplicateKey(String::from(
                key::SAMPLES_WITH_DATA_COUNT
            )))
        );

        Ok(())
    }
}
