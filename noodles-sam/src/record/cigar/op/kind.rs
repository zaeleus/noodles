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
    SequenceMatch,
    /// A sequence mismatch (`X`).
    SequenceMismatch,
}

impl Kind {
    /// Returns whether the operation kind causes the alignment to consume the read.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::cigar::op::Kind;
    /// assert!(Kind::Match.consumes_read());
    /// assert!(Kind::Insertion.consumes_read());
    /// assert!(!Kind::Deletion.consumes_read());
    /// ```
    pub fn consumes_read(&self) -> bool {
        matches!(
            self,
            Self::Match
                | Self::Insertion
                | Self::SoftClip
                | Self::SequenceMatch
                | Self::SequenceMismatch
        )
    }

    /// Returns whether the operation kind causes the alignment to consume the reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::cigar::op::Kind;
    /// assert!(Kind::Match.consumes_reference());
    /// assert!(!Kind::Insertion.consumes_reference());
    /// assert!(Kind::Deletion.consumes_reference());
    /// ```
    pub fn consumes_reference(&self) -> bool {
        matches!(
            self,
            Self::Match
                | Self::Deletion
                | Self::Skip
                | Self::SequenceMatch
                | Self::SequenceMismatch
        )
    }
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
            "=" => Ok(Self::SequenceMatch),
            "X" => Ok(Self::SequenceMismatch),
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
            Kind::SequenceMatch => '=',
            Kind::SequenceMismatch => 'X',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_consumes_read() {
        assert!(Kind::Match.consumes_read());
        assert!(Kind::Insertion.consumes_read());
        assert!(!Kind::Deletion.consumes_read());
        assert!(!Kind::Skip.consumes_read());
        assert!(Kind::SoftClip.consumes_read());
        assert!(!Kind::HardClip.consumes_read());
        assert!(!Kind::Pad.consumes_read());
        assert!(Kind::SequenceMatch.consumes_read());
        assert!(Kind::SequenceMismatch.consumes_read());
    }

    #[test]
    fn test_consumes_reference() {
        assert!(Kind::Match.consumes_reference());
        assert!(!Kind::Insertion.consumes_reference());
        assert!(Kind::Deletion.consumes_reference());
        assert!(Kind::Skip.consumes_reference());
        assert!(!Kind::SoftClip.consumes_reference());
        assert!(!Kind::HardClip.consumes_reference());
        assert!(!Kind::Pad.consumes_reference());
        assert!(Kind::SequenceMatch.consumes_reference());
        assert!(Kind::SequenceMismatch.consumes_reference());
    }

    #[test]
    fn test_fmt() {
        assert_eq!(Kind::Match.to_string(), "M");
        assert_eq!(Kind::Insertion.to_string(), "I");
        assert_eq!(Kind::Deletion.to_string(), "D");
        assert_eq!(Kind::Skip.to_string(), "N");
        assert_eq!(Kind::SoftClip.to_string(), "S");
        assert_eq!(Kind::HardClip.to_string(), "H");
        assert_eq!(Kind::Pad.to_string(), "P");
        assert_eq!(Kind::SequenceMatch.to_string(), "=");
        assert_eq!(Kind::SequenceMismatch.to_string(), "X");
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
        assert_eq!("=".parse(), Ok(Kind::SequenceMatch));
        assert_eq!("X".parse(), Ok(Kind::SequenceMismatch));

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
        assert_eq!(char::from(Kind::SequenceMatch), '=');
        assert_eq!(char::from(Kind::SequenceMismatch), 'X');
    }
}
