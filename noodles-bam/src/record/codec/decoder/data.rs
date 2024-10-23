pub mod field;

pub(crate) use self::field::read_field;

use std::{error, fmt};

use noodles_sam::alignment::{record::data::field::Tag, record_buf::Data};

/// An error when raw BAM record data fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// A tag is duplicated.
    DuplicateTag(Tag),
    /// A field is invalid.
    InvalidField(field::DecodeError),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::DuplicateTag(_) => None,
            Self::InvalidField(e) => Some(e),
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag:?}"),
            Self::InvalidField(e) => {
                write!(f, "invalid field")?;

                if let Some(tag) = e.tag() {
                    write!(f, ": {tag:?}")?;
                }

                Ok(())
            }
        }
    }
}

pub(crate) fn read_data(src: &mut &[u8], data: &mut Data) -> Result<(), DecodeError> {
    data.clear();

    while !src.is_empty() {
        let (tag, value) = read_field(src).map_err(DecodeError::InvalidField)?;

        if let Some((t, _)) = data.insert(tag, value) {
            return Err(DecodeError::DuplicateTag(t));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_data() -> Result<(), DecodeError> {
        use noodles_sam::alignment::record_buf::data::field::Value;

        fn t(mut src: &[u8], actual: &mut Data, expected: &Data) -> Result<(), DecodeError> {
            read_data(&mut src, actual)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        let mut buf = Data::default();

        let expected = Data::default();
        t(&[], &mut buf, &expected)?;

        let expected = [(Tag::ALIGNMENT_HIT_COUNT, Value::UInt8(1))]
            .into_iter()
            .collect();

        t(
            &[b'N', b'H', b'C', 0x01], // NH:C:1
            &mut buf,
            &expected,
        )?;

        let expected = [
            (Tag::ALIGNMENT_HIT_COUNT, Value::UInt8(1)),
            (Tag::READ_GROUP, Value::from("rg0")),
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
            read_data(&mut src, &mut buf),
            Err(DecodeError::DuplicateTag(Tag::ALIGNMENT_HIT_COUNT))
        );

        Ok(())
    }
}
