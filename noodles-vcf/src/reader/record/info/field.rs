pub(super) mod value;

use self::value::parse_value;
use super::ParseError;
use crate::{
    header::record::value::map::info::Type,
    reader::record::MISSING,
    record::info::field::{Key, Value},
    Header,
};

pub(super) fn parse_field(header: &Header, s: &str) -> Result<(Key, Option<Value>), ParseError> {
    use crate::header::record::value::map::info::definition::definition;

    const MAX_COMPONENTS: usize = 2;
    const SEPARATOR: char = '=';

    let mut components = s.splitn(MAX_COMPONENTS, SEPARATOR);

    let raw_key = components.next().unwrap_or_default();
    let key = raw_key.parse().map_err(ParseError::InvalidKey)?;

    let (number, ty) = header
        .infos()
        .get(&key)
        .map(|info| (info.number(), info.ty()))
        .or_else(|| definition(header.file_format(), &key).map(|(n, t, _)| (n, t)))
        .unwrap_or_default();

    let raw_value = components.next();

    let value = if matches!(ty, Type::Flag) {
        match raw_value.unwrap_or_default() {
            MISSING => None,
            t => parse_value(number, ty, t)
                .map(Some)
                .map_err(ParseError::InvalidValue)?,
        }
    } else if matches!(key, Key::Other(_)) {
        match raw_value {
            Some(MISSING) => None,
            Some(t) => parse_value(number, ty, t)
                .map(Some)
                .map_err(ParseError::InvalidValue)?,
            None => Some(Value::Flag),
        }
    } else if let Some(t) = raw_value {
        match t {
            MISSING => None,
            _ => parse_value(number, ty, t)
                .map(Some)
                .map_err(ParseError::InvalidValue)?,
        }
    } else {
        return Err(ParseError::MissingValue);
    };

    Ok((key, value))
}
