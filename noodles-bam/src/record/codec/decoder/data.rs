pub mod field;

pub(crate) use self::field::get_field;

use std::{error, fmt};

use bytes::Buf;
use noodles_sam::record::{data::field::Tag, Data};

/// An error when raw BAM record data fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A tag is duplicated.
    DuplicateTag(Tag),
    /// A field is invalid.
    InvalidField(field::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::DuplicateTag(_) => None,
            Self::InvalidField(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
            Self::InvalidField(e) => {
                write!(f, "invalid field")?;

                if let Some(tag) = e.tag() {
                    write!(f, ": {tag}")?;
                }

                Ok(())
            }
        }
    }
}

pub(crate) fn get_data<B>(src: &mut B, data: &mut Data) -> Result<(), ParseError>
where
    B: Buf,
{
    data.clear();

    while src.has_remaining() {
        let (tag, value) = get_field(src).map_err(ParseError::InvalidField)?;

        if let Some((t, _)) = data.insert(tag, value) {
            return Err(ParseError::DuplicateTag(t));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_get_data() -> Result<(), ParseError> {
        use noodles_sam::record::data::field::{tag, Value};

        fn t(mut src: &[u8], actual: &mut Data, expected: &Data) -> Result<(), ParseError> {
            get_data(&mut src, actual)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        let mut buf = Data::default();

        let expected = Data::default();
        t(&[], &mut buf, &expected)?;

        let expected = [(tag::ALIGNMENT_HIT_COUNT, Value::UInt8(1))]
            .into_iter()
            .collect();

        t(
            &[b'N', b'H', b'C', 0x01], // NH:C:1
            &mut buf,
            &expected,
        )?;

        let expected = [
            (tag::ALIGNMENT_HIT_COUNT, Value::UInt8(1)),
            (tag::READ_GROUP, Value::String(String::from("rg0"))),
        ]
        .into_iter()
        .collect();

        t(
            &[
                b'N', b'H', b'C', 0x01, // NH:C:1
                b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
            ],
            &mut buf,
            &expected,
        )?;

        let data = [
            b'N', b'H', b'C', 0x01, // NH:C:1
            b'N', b'H', b'C', 0x01, // NH:C:1
        ];
        let mut src = &data[..];
        assert_eq!(
            get_data(&mut src, &mut buf),
            Err(ParseError::DuplicateTag(tag::ALIGNMENT_HIT_COUNT))
        );

        Ok(())
    }
}
