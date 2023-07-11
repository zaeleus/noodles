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
    if let Some(modifications) = parse_short_codes(src)? {
        Ok(modifications)
    } else if let Some(modification) = parse_chebi_id(src)? {
        Ok(vec![modification])
    } else {
        Err(ParseError::Invalid)
    }
}

fn parse_short_codes(src: &mut &[u8]) -> Result<Option<Vec<Modification>>, ParseError> {
    let raw_codes = take_while(src, |b| b.is_ascii_lowercase());

    if raw_codes.is_empty() {
        return Ok(None);
    }

    let modifications = raw_codes
        .iter()
        .copied()
        .map(|raw_code| Modification::try_from(raw_code).map_err(|_| ParseError::Invalid))
        .collect::<Result<_, _>>()?;

    Ok(Some(modifications))
}

fn parse_chebi_id(src: &mut &[u8]) -> Result<Option<Modification>, ParseError> {
    let raw_id = take_while(src, |b| b.is_ascii_digit());

    if raw_id.is_empty() {
        Ok(None)
    } else {
        lexical_core::parse(raw_id)
            .map(Modification::ChebiId)
            .map(Some)
            .map_err(|_| ParseError::Invalid)
    }
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
        use crate::record::data::field::value::base_modifications::group::modification;

        let mut src = &b"m"[..];
        assert_eq!(
            parse_modifications(&mut src),
            Ok(vec![modification::FIVE_METHYLCYTOSINE])
        );

        let mut src = &b"mh"[..];
        assert_eq!(
            parse_modifications(&mut src),
            Ok(vec![
                modification::FIVE_METHYLCYTOSINE,
                modification::FIVE_HYDROXYMETHYLCYTOSINE,
            ])
        );

        let mut src = &b"27551"[..];
        assert_eq!(
            parse_modifications(&mut src),
            Ok(vec![Modification::ChebiId(27551)])
        );

        let mut src = &b""[..];
        assert_eq!(parse_modifications(&mut src), Err(ParseError::Invalid));
    }
}
