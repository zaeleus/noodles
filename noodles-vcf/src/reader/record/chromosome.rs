use std::{error, fmt};

use noodles_core as core;

use crate::record::Chromosome;

/// An error when a raw VCF record chromosome fails to parse.
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

impl From<ParseError> for core::Error {
    fn from(e: ParseError) -> Self {
        Self::new(core::error::Kind::Parse, e)
    }
}

pub(super) fn parse_chromosome(s: &str, chromosome: &mut Chromosome) -> Result<(), ParseError> {
    // symbol
    if let Some(t) = s.strip_prefix('<') {
        if let Some(t) = t.strip_suffix('>') {
            if !matches!(chromosome, Chromosome::Symbol(symbol) if symbol == t) {
                *chromosome = Chromosome::Symbol(t.into());
            }

            return Ok(());
        }
    }

    // name
    if !matches!(chromosome, Chromosome::Name(name) if name == s) {
        if is_valid_name(s) {
            *chromosome = Chromosome::Name(s.into());
        } else {
            return Err(ParseError::Invalid);
        }
    }

    Ok(())
}

// § 1.4.7 "Contig field format"
fn is_valid_name_char(c: char) -> bool {
    ('!'..='~').contains(&c)
        && !matches!(
            c,
            '\\' | ',' | '"' | '`' | '\'' | '(' | ')' | '[' | ']' | '{' | '}' | '<' | '>',
        )
}

fn is_valid_name(s: &str) -> bool {
    let mut chars = s.chars();

    let is_valid_first_char = chars
        .next()
        .map(|c| c != '*' && c != '=' && is_valid_name_char(c))
        .unwrap_or_default();

    is_valid_first_char && chars.all(is_valid_name_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_chromosome() -> Result<(), ParseError> {
        let mut chromosome = Chromosome::Name(String::from("."));

        parse_chromosome("sq0", &mut chromosome)?;
        assert_eq!(chromosome, Chromosome::Name(String::from("sq0")));

        parse_chromosome("<sq0>", &mut chromosome)?;
        assert_eq!(chromosome, Chromosome::Symbol(String::from("sq0")));

        assert_eq!(
            parse_chromosome("", &mut chromosome),
            Err(ParseError::Invalid)
        );

        Ok(())
    }

    #[test]
    fn test_is_valid_name() {
        assert!(is_valid_name("sq0"));

        assert!(!is_valid_name(""));
        assert!(!is_valid_name("sq 0"));
        assert!(!is_valid_name("sq[0]"));
        assert!(!is_valid_name(">sq0"));
        assert!(!is_valid_name("*sq0"));
        assert!(!is_valid_name("=sq0"));
    }
}
