mod kind;

use std::str::FromStr;

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

impl FromStr for Record {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut pieces = s.split(DELIMITER);

        let kind = pieces
            .next()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| ())?;

        let record = if let Kind::Comment = kind {
            let comment = pieces.next().unwrap();
            Self::Comment(comment.into())
        } else {
            let fields = pieces
                .map(|field| {
                    let mut field_pieces = field.splitn(2, DATA_FIELD_DELIMITER);
                    let tag = field_pieces.next().ok_or_else(|| ())?;
                    let value = field_pieces.next().ok_or_else(|| ())?;
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
