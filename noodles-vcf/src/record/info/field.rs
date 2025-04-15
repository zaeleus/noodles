mod value;

use std::io;

use self::value::parse_value;
use crate::{header::record::value::map::info::Type, variant::record::info::field::Value, Header};

pub(super) fn parse_field<'a>(
    buf: &'a str,
    header: &Header,
) -> io::Result<(&'a str, Option<Value<'a>>)> {
    use crate::header::record::value::map::info::definition::definition;

    const MAX_COMPONENTS: usize = 2;
    const MISSING: &str = ".";
    const SEPARATOR: char = '=';

    if buf.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "empty info field",
        ));
    }

    let mut components = buf.splitn(MAX_COMPONENTS, SEPARATOR);
    let key = components.next().unwrap_or_default();
    let raw_value = components.next();

    let definition = header
        .infos()
        .get(key)
        .map(|info| (info.number(), info.ty()))
        .or_else(|| definition(header.file_format(), key).map(|(n, t, _)| (n, t)));

    if definition.is_none() && raw_value.is_none() {
        return Ok((key, Some(Value::Flag)));
    }

    let (number, ty) = definition.unwrap_or_default();

    let value = match raw_value {
        Some(MISSING) => None,
        Some(t) => parse_value(t, number, ty).map(Some)?,
        None if ty == Type::Flag => Some(Value::Flag),
        None => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "missing info field value",
            ))
        }
    };

    Ok((key, value))
}
