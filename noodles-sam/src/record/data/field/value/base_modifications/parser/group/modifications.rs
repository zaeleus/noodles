use std::{error, fmt};

use crate::record::data::field::value::base_modifications::group::Modification;

/// An error returned when a base modifications group modifications fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

pub fn parse_modifications(src: &mut &[u8]) -> Result<Vec<Modification>, ParseError> {
    parse_short_codes(src)
}

fn parse_short_codes(src: &mut &[u8]) -> Result<Vec<Modification>, ParseError> {
    let raw_codes = take_while(src, |b| b.is_ascii_lowercase());

    if raw_codes.is_empty() {
        return Err(ParseError::Invalid);
    }

    let codes = raw_codes
        .iter()
        .copied()
        .map(|raw_code| Modification::try_from(raw_code).map_err(|_| ParseError::Invalid))
        .collect::<Result<_, _>>()?;

    Ok(codes)
}

fn take_while<'a, P>(src: &mut &'a [u8], predicate: P) -> &'a [u8]
where
    P: Fn(&u8) -> bool,
{
    match src.iter().position(|b| !predicate(b)) {
        Some(i) => {
            let (token, rest) = src.split_at(i);
            *src = rest;
            token
        }
        None => src,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_modifications() {
        let mut src = &b"m"[..];
        assert_eq!(
            parse_modifications(&mut src),
            Ok(vec![Modification::FiveMethylcytosine])
        );

        let mut src = &b"mh"[..];
        assert_eq!(
            parse_modifications(&mut src),
            Ok(vec![
                Modification::FiveMethylcytosine,
                Modification::FiveHydroxymethylcytosine
            ])
        );

        // TODO
        let mut src = &b"27551"[..];
        assert_eq!(parse_modifications(&mut src), Err(ParseError::Invalid));

        let mut src = &b""[..];
        assert_eq!(parse_modifications(&mut src), Err(ParseError::Invalid));
    }
}
