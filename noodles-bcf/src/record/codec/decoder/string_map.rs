use std::{error, fmt, num};

use crate::{
    header::string_maps::StringMap,
    record::codec::{
        value::{Array, Int16, Int32, Int8},
        Value,
    },
};

use super::value::read_value;

pub fn read_string_map_index(src: &mut &[u8]) -> Result<usize, DecodeError> {
    let value = read_value(src).map_err(DecodeError::InvalidValue)?;

    match value.as_ref().and_then(|v| v.as_int()) {
        Some(i) => usize::try_from(i).map_err(DecodeError::InvalidIndex),
        None => Err(DecodeError::InvalidIndexValue),
    }
}

pub fn read_string_map_indices(src: &mut &[u8]) -> Result<Vec<usize>, DecodeError> {
    let value = read_value(src).map_err(DecodeError::InvalidValue)?;

    let indices = match value {
        Some(Value::Int8(Some(Int8::Value(i)))) => vec![i32::from(i)],
        Some(Value::Array(Array::Int8(values))) => values
            .iter()
            .map(|result| result.map(i32::from))
            .collect::<Result<_, _>>()
            .map_err(|_| DecodeError::UnexpectedEof)?,
        Some(Value::Int16(Some(Int16::Value(i)))) => vec![i32::from(i)],
        Some(Value::Array(Array::Int16(values))) => values
            .iter()
            .map(|result| result.map(i32::from))
            .collect::<Result<_, _>>()
            .map_err(|_| DecodeError::UnexpectedEof)?,
        Some(Value::Int32(Some(Int32::Value(i)))) => vec![i],
        Some(Value::Array(Array::Int32(values))) => values
            .iter()
            .collect::<Result<_, _>>()
            .map_err(|_| DecodeError::UnexpectedEof)?,
        None => Vec::new(),
        _ => return Err(DecodeError::InvalidIndexValue),
    };

    indices
        .into_iter()
        .map(|i| usize::try_from(i).map_err(DecodeError::InvalidIndex))
        .collect()
}

pub fn read_string_map_entry<'m>(
    src: &mut &[u8],
    string_map: &'m StringMap,
) -> Result<&'m str, DecodeError> {
    read_string_map_index(src)
        .and_then(|j| string_map.get_index(j).ok_or(DecodeError::MissingEntry))
}

#[allow(clippy::enum_variant_names)]
#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    UnexpectedEof,
    InvalidValue(super::value::DecodeError),
    InvalidIndex(num::TryFromIntError),
    InvalidIndexValue,
    MissingEntry,
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidValue(e) => Some(e),
            Self::InvalidIndex(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidValue(_) => write!(f, "invalid value"),
            Self::InvalidIndex(_) => write!(f, "invalid index"),
            Self::InvalidIndexValue => write!(f, "invalid index value"),
            Self::MissingEntry => write!(f, "missing entry"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_string_map_index() {
        fn t(mut src: &[u8], expected: usize) {
            assert_eq!(read_string_map_index(&mut src), Ok(expected));
        }

        // Some(Type::Int8(Some(Int8::Value(8))))
        t(&[0x11, 0x08], 8);

        // Some(Type::Int16(Some(Int16::Value(13))))
        t(&[0x12, 0x0d, 0x00], 13);

        // Some(Type::Int32(Some(Int32::Value(21))))
        t(&[0x13, 0x15, 0x00, 0x00, 0x00], 21);

        // Some(Type::Int8(Some(Int8::Value(-5))))
        let mut src = &[0x11, 0xfb][..];
        assert!(matches!(
            read_string_map_index(&mut src),
            Err(DecodeError::InvalidIndex(_))
        ));

        // Some(Type::String(Some(String::from("n"))))
        let mut src = &[0x17, b'n'][..];
        assert_eq!(
            read_string_map_index(&mut src),
            Err(DecodeError::InvalidIndexValue)
        );
    }

    #[test]
    fn test_read_string_map_indices() {
        fn t(mut src: &[u8], expected: &[usize]) {
            assert_eq!(read_string_map_indices(&mut src), Ok(expected.into()));
        }

        // None
        t(&[0x00], &[]);

        // Some(Type::Int8(Some(Int8::Value(2))))
        t(&[0x11, 0x02], &[2]);
        // Some(Type::Array(Array::Int8(vec![2, 3])))
        t(&[0x21, 0x02, 0x03], &[2, 3]);

        // Some(Type::Int16(Some(Int16::Value(5))))
        t(&[0x12, 0x05, 0x00], &[5]);
        // Some(Type::Array(Array::Int16(vec![5, 8])))
        t(&[0x22, 0x05, 0x00, 0x08, 0x00], &[5, 8]);

        // Some(Type::Int32(Some(Int32::Value(13))))
        t(&[0x13, 0x0d, 0x00, 0x00, 0x00], &[13]);
        // Some(Type::Array(Array::Int32(vec![13, 21])))
        t(
            &[0x23, 0x0d, 0x00, 0x00, 0x00, 0x15, 0x00, 0x00, 0x00],
            &[13, 21],
        );

        // Some(Type::Int8(Some(Int8::Value(-5))))
        let mut src = &[0x11, 0xfb][..];
        assert!(matches!(
            read_string_map_indices(&mut src),
            Err(DecodeError::InvalidIndex(_))
        ));

        // Some(Type::Int8(None))
        let mut src = &[0x01][..];
        assert_eq!(
            read_string_map_indices(&mut src),
            Err(DecodeError::InvalidIndexValue)
        );

        // Some(Type::Int16(None))
        let mut src = &[0x02][..];
        assert_eq!(
            read_string_map_indices(&mut src),
            Err(DecodeError::InvalidIndexValue)
        );

        // Some(Type::Int32(None))
        let mut src = &[0x03][..];
        assert_eq!(
            read_string_map_indices(&mut src),
            Err(DecodeError::InvalidIndexValue)
        );

        // Some(Type::String(Some(String::from("n"))))
        let mut src = &[0x17, b'n'][..];
        assert_eq!(
            read_string_map_indices(&mut src),
            Err(DecodeError::InvalidIndexValue)
        );
    }

    #[test]
    fn test_read_string_map_entry() {
        let mut string_map = StringMap::default();
        string_map.insert(String::from("PASS"));
        string_map.insert(String::from("NS"));
        string_map.insert(String::from("DP"));

        // Some(Type::Int8(Some(Int8::Value(1))))
        let mut src = &[0x11, 0x01][..];
        assert_eq!(read_string_map_entry(&mut src, &string_map), Ok("NS"));

        // Some(Type::Int8(Some(Int8::Value(8))))
        let mut src = &[0x11, 0x08][..];
        assert_eq!(
            read_string_map_entry(&mut src, &string_map),
            Err(DecodeError::MissingEntry)
        );
    }
}
