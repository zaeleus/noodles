pub(crate) mod header;
mod tag;

use std::str;

pub(crate) use self::header::parse_header;
use self::tag::parse_tag;

fn parse_string<'a>(src: &mut &'a [u8]) -> Result<&'a str, str::Utf8Error> {
    use memchr::memchr;

    const DELIMITER: u8 = b'\t';

    let i = memchr(DELIMITER, src).unwrap_or(src.len());
    let (buf, rest) = src.split_at(i);

    *src = rest;

    str::from_utf8(buf)
}
