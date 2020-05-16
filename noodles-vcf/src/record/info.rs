pub mod field;

pub use self::field::Field;

use std::{error, fmt, ops::Deref, str::FromStr};

use super::MISSING_FIELD;

const DELIMITER: char = ';';

#[derive(Debug, Default)]
pub struct Info(Vec<Field>);

impl Deref for Info {
    type Target = [Field];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Info {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_empty() {
            f.write_str(MISSING_FIELD)
        } else {
            for (i, field) in self.iter().enumerate() {
                if i > 0 {
                    write!(f, "{}", DELIMITER)?
                }

                write!(f, "{}", field)?;
            }

            Ok(())
        }
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid info: {}", self.0)
    }
}

impl FromStr for Info {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == MISSING_FIELD {
            return Ok(Info::default());
        }

        s.split(DELIMITER)
            .map(|s| s.parse())
            .collect::<Result<_, _>>()
            .map(Info)
            .map_err(|_| ParseError(s.into()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let info = Info::default();
        assert_eq!(info.to_string(), ".");

        let info = Info(vec![Field::new(
            field::Key::SamplesWithDataCount,
            field::Value::Integer(2),
        )]);
        assert_eq!(info.to_string(), "NS=2");

        let info = Info(vec![
            Field::new(field::Key::SamplesWithDataCount, field::Value::Integer(2)),
            Field::new(
                field::Key::AlleleFrequencies,
                field::Value::FloatArray(vec![0.333, 0.667]),
            ),
        ]);
        assert_eq!(info.to_string(), "NS=2;AF=0.333,0.667");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let actual: Info = ".".parse()?;
        assert!(actual.is_empty());

        let actual: Info = "NS=2".parse()?;
        assert_eq!(actual.len(), 1);

        let actual: Info = "NS=2;AF=0.333,0.667".parse()?;
        assert_eq!(actual.len(), 2);

        Ok(())
    }
}
