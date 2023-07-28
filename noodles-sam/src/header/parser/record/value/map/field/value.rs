use std::str;

pub fn parse_value<'a>(src: &mut &'a [u8]) -> Result<&'a str, str::Utf8Error> {
    use memchr::memchr;

    const DELIMITER: u8 = b'\t';

    let i = memchr(DELIMITER, src).unwrap_or(src.len());
    let (buf, rest) = src.split_at(i);

    *src = rest;

    str::from_utf8(buf)
}
