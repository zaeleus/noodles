mod kind;

use std::{error, fmt, str::FromStr};

pub use self::kind::Kind;

type Field = (String, String);

const DELIMITER: char = '\t';
const DATA_FIELD_DELIMITER: char = ':';

#[derive(Debug)]
pub enum Record {
    Header(Vec<Field>),
    ReferenceSequence(Vec<Field>),
    ReadGroup(Vec<Field>),
    Program(Vec<Field>),
    Comment(String),
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid record: {}", self.0)
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut pieces = s.split(DELIMITER);

        let kind = pieces
            .next()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| ParseError(s.into()))?;

        let record = if let Kind::Comment = kind {
            let comment = pieces.next().unwrap();
            Self::Comment(comment.into())
        } else {
            let fields = pieces
                .map(|field| {
                    let mut field_pieces = field.splitn(2, DATA_FIELD_DELIMITER);
                    let tag = field_pieces.next().ok_or_else(|| ParseError(s.into()))?;
                    let value = field_pieces.next().ok_or_else(|| ParseError(s.into()))?;
                    Ok((tag.into(), value.into()))
                })
                .collect::<Result<_, _>>()?;

            match kind {
                Kind::Header => Self::Header(fields),
                Kind::ReferenceSequence => Self::ReferenceSequence(fields),
                Kind::ReadGroup => Self::ReadGroup(fields),
                Kind::Program => Self::Program(fields),
                Kind::Comment => unreachable!(),
            }
        };

        Ok(record)
    }
}
