use super::ParseError;
use crate::header::record::value::map::header::Version;

pub(super) fn parse_version(src: &[u8]) -> Result<Version, ParseError> {
    const DELIMITER: u8 = b'.';

    fn split_once(buf: &[u8], delimiter: u8) -> Option<(&[u8], &[u8])> {
        let i = buf.iter().position(|&b| b == delimiter)?;
        Some((&buf[..i], &buf[i + 1..]))
    }

    match split_once(src, DELIMITER) {
        Some((a, b)) => {
            let major = lexical_core::parse(a).map_err(|_| ParseError::InvalidVersion)?;
            let minor = lexical_core::parse(b).map_err(|_| ParseError::InvalidVersion)?;
            Ok(Version::new(major, minor))
        }
        None => Err(ParseError::InvalidVersion),
    }
}
