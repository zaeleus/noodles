//! Base modifications group modification.

use std::{error, fmt};

/// 5-Methylcytosine (m; 5mC; CHEBI:27551).
pub const FIVE_METHYLCYTOSINE: Modification = Modification::Code(b'm');

/// 5-Hydroxymethylcytosine (h; 5hmC; CHEBI:76792).
pub const FIVE_HYDROXYMETHYLCYTOSINE: Modification = Modification::Code(b'h');

/// 5-Formylcytosine (f; 5fC; CHEBI:76794).
pub const FIVE_FORMYLCYTOSINE: Modification = Modification::Code(b'f');

/// 5-Carboxylcytosine (c; 5caC; CHEBI:76793).
pub const FIVE_CARBOXYLCYTOSINE: Modification = Modification::Code(b'c');

/// 5-Hydroxymethyluracil (g; 5hmU; CHEBI:16964).
pub const FIVE_HYDROXYMETHYLURACIL: Modification = Modification::Code(b'g');

/// 5-Formyluracil (e; 5fU; CHEBI:80961).
pub const FIVE_FORMYLURACIL: Modification = Modification::Code(b'e');

/// 5-Carboxyluracil (b; 5caU; CHEBI:17477).
pub const FIVE_CARBOXYLURACIL: Modification = Modification::Code(b'b');

/// 6-Methyladenine (a; 6mA; CHEBI:28871).
pub const SIX_METHYLADENINE: Modification = Modification::Code(b'a');

/// 8-Oxoguanine (o; 8oxoG; CHEBI:44605).
pub const EIGHT_OXOGUANINE: Modification = Modification::Code(b'o');

/// Xanthosine (n; Xao; CHEBI:18107).
pub const XANTHOSINE: Modification = Modification::Code(b'n');

/// A base modification.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Modification {
    /// A modification code.
    Code(u8),
    /// A ChEBI ID.
    ChebiId(u32),
}

/// An error returned when a base modifications group modification fails to parse.
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

impl TryFrom<u8> for Modification {
    type Error = ParseError;

    fn try_from(b: u8) -> Result<Self, Self::Error> {
        if b.is_ascii_lowercase() || is_ambiguity_code(b) {
            Ok(Self::Code(b))
        } else {
            Err(ParseError::Invalid)
        }
    }
}

fn is_ambiguity_code(b: u8) -> bool {
    matches!(b, b'C' | b'T' | b'U' | b'A' | b'G' | b'N')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_u8_for_modification() {
        fn t(b: u8, expected: Modification) {
            assert_eq!(Modification::try_from(b), Ok(expected));
        }

        t(b'm', FIVE_METHYLCYTOSINE);
        t(b'h', FIVE_HYDROXYMETHYLCYTOSINE);
        t(b'f', FIVE_FORMYLCYTOSINE);
        t(b'c', FIVE_CARBOXYLCYTOSINE);
        t(b'C', Modification::Code(b'C'));

        t(b'g', FIVE_HYDROXYMETHYLURACIL);
        t(b'e', FIVE_FORMYLURACIL);
        t(b'b', FIVE_CARBOXYLURACIL);
        t(b'T', Modification::Code(b'T'));

        t(b'U', Modification::Code(b'U'));

        t(b'a', SIX_METHYLADENINE);
        t(b'A', Modification::Code(b'A'));

        t(b'o', EIGHT_OXOGUANINE);
        t(b'G', Modification::Code(b'G'));

        t(b'n', XANTHOSINE);
        t(b'N', Modification::Code(b'N'));

        t(b'z', Modification::Code(b'z'));

        assert_eq!(Modification::try_from(b'?'), Err(ParseError::Invalid));
    }
}
