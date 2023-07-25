pub(crate) mod header;
pub(crate) mod reference_sequence;
mod tag;

pub(crate) use self::{header::parse_header, reference_sequence::parse_reference_sequence};

use std::str;

use self::tag::parse_tag;

fn parse_string<'a>(src: &mut &'a [u8]) -> Result<&'a str, str::Utf8Error> {
    use memchr::memchr;

    const DELIMITER: u8 = b'\t';

    let i = memchr(DELIMITER, src).unwrap_or(src.len());
    let (buf, rest) = src.split_at(i);

    *src = rest;

    str::from_utf8(buf)
}
