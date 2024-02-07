mod value;

use std::io;

use self::value::parse_value;
use crate::{header::record::value::map::info::Type, variant::record::info::field::Value, Header};

pub(super) fn parse_field<'a>(
    src: &mut &'a str,
    header: &Header,
) -> io::Result<(&'a str, Option<Value<'a>>)> {
    use crate::header::record::value::map::info::definition::definition;

    const DELIMITER: char = ';';
    const MAX_COMPONENTS: usize = 2;
    const MISSING: &str = ".";
    const SEPARATOR: char = '=';

    let (buf, rest) = match src.find(DELIMITER) {
        Some(i) => {
            let (buf, rest) = src.split_at(i);
            (buf, &rest[1..])
        }
        None => src.split_at(src.len()),
    };

    *src = rest;

    if buf.is_empty() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
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
        None => return Err(io::Error::new(io::ErrorKind::InvalidData, "missing value")),
    };

    Ok((key, value))
}
