//! SAM record CIGAR operation kind.

use std::{
    error,
    fmt::{self, Write},
    str::FromStr,
};

/// A SAM record CIGAR operation kind.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    /// An alignment match (`M`).
    Match,
    /// An insertion into the reference (`I`).
    Insertion,
    /// A deletion from the reference (`D`).
    Deletion,
    /// A skipped region from the reference (`N`).
    Skip,
    /// A soft clip (`S`).
    SoftClip,
    /// A hard clip (`H`).
    HardClip,
    /// Padding (`P`).
    Pad,
    /// A sequence match (`=`).
    SeqMatch,
    /// A sequence mismatch (`X`).
    SeqMismatch,
}

impl fmt::Display for Kind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_char(char::from(*self))
    }
}

/// An error returned when a raw SAM record CIGAR operation kind fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl FromStr for Kind {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "M" => Ok(Self::Match),
            "I" => Ok(Self::Insertion),
            "D" => Ok(Self::Deletion),
            "N" => Ok(Self::Skip),
            "S" => Ok(Self::SoftClip),
            "H" => Ok(Self::HardClip),
            "P" => Ok(Self::Pad),
            "=" => Ok(Self::SeqMatch),
            "X" => Ok(Self::SeqMismatch),
            _ => Err(ParseError(s.into())),
        }
    }
}

impl From<Kind> for char {
    fn from(kind: Kind) -> Self {
        match kind {
            Kind::Match => 'M',
            Kind::Insertion => 'I',
            Kind::Deletion => 'D',
            Kind::Skip => 'N',
            Kind::SoftClip => 'S',
            Kind::HardClip => 'H',
            Kind::Pad => 'P',
            Kind::SeqMatch => '=',
            Kind::SeqMismatch => 'X',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Kind::Match.to_string(), "M");
        assert_eq!(Kind::Insertion.to_string(), "I");
        assert_eq!(Kind::Deletion.to_string(), "D");
        assert_eq!(Kind::Skip.to_string(), "N");
        assert_eq!(Kind::SoftClip.to_string(), "S");
        assert_eq!(Kind::HardClip.to_string(), "H");
        assert_eq!(Kind::Pad.to_string(), "P");
        assert_eq!(Kind::SeqMatch.to_string(), "=");
        assert_eq!(Kind::SeqMismatch.to_string(), "X");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("M".parse(), Ok(Kind::Match));
        assert_eq!("I".parse(), Ok(Kind::Insertion));
        assert_eq!("D".parse(), Ok(Kind::Deletion));
        assert_eq!("N".parse(), Ok(Kind::Skip));
        assert_eq!("S".parse(), Ok(Kind::SoftClip));
        assert_eq!("H".parse(), Ok(Kind::HardClip));
        assert_eq!("P".parse(), Ok(Kind::Pad));
        assert_eq!("=".parse(), Ok(Kind::SeqMatch));
        assert_eq!("X".parse(), Ok(Kind::SeqMismatch));

        assert_eq!("".parse::<Kind>(), Err(ParseError(String::from(""))));
        assert_eq!("O".parse::<Kind>(), Err(ParseError(String::from("O"))));
    }

    #[test]
    fn test_from_kind_for_char() {
        assert_eq!(char::from(Kind::Match), 'M');
        assert_eq!(char::from(Kind::Insertion), 'I');
        assert_eq!(char::from(Kind::Deletion), 'D');
        assert_eq!(char::from(Kind::Skip), 'N');
        assert_eq!(char::from(Kind::SoftClip), 'S');
        assert_eq!(char::from(Kind::HardClip), 'H');
        assert_eq!(char::from(Kind::Pad), 'P');
        assert_eq!(char::from(Kind::SeqMatch), '=');
        assert_eq!(char::from(Kind::SeqMismatch), 'X');
    }
}
