use std::{error, fmt, str::FromStr};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    Match,
    Insertion,
    Deletion,
    Skip,
    SoftClip,
    HardClip,
    Pad,
    SeqMatch,
    SeqMismatch,
}

impl Kind {
    pub fn symbol(self) -> char {
        match self {
            Self::Match => 'M',
            Self::Insertion => 'I',
            Self::Deletion => 'D',
            Self::Skip => 'N',
            Self::SoftClip => 'S',
            Self::HardClip => 'H',
            Self::Pad => 'P',
            Self::SeqMatch => '=',
            Self::SeqMismatch => 'X',
        }
    }
}

impl fmt::Display for Kind {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid cigar op kind: {}", self.0)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_symbol() {
        assert_eq!(Kind::Match.symbol(), 'M');
        assert_eq!(Kind::Insertion.symbol(), 'I');
        assert_eq!(Kind::Deletion.symbol(), 'D');
        assert_eq!(Kind::Skip.symbol(), 'N');
        assert_eq!(Kind::SoftClip.symbol(), 'S');
        assert_eq!(Kind::HardClip.symbol(), 'H');
        assert_eq!(Kind::Pad.symbol(), 'P');
        assert_eq!(Kind::SeqMatch.symbol(), '=');
        assert_eq!(Kind::SeqMismatch.symbol(), 'X');
    }

    #[test]
    fn test_fmt() {
        assert_eq!(format!("{}", Kind::Match), "M");
        assert_eq!(format!("{}", Kind::Insertion), "I");
        assert_eq!(format!("{}", Kind::Deletion), "D");
        assert_eq!(format!("{}", Kind::Skip), "N");
        assert_eq!(format!("{}", Kind::SoftClip), "S");
        assert_eq!(format!("{}", Kind::HardClip), "H");
        assert_eq!(format!("{}", Kind::Pad), "P");
        assert_eq!(format!("{}", Kind::SeqMatch), "=");
        assert_eq!(format!("{}", Kind::SeqMismatch), "X");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("M".parse::<Kind>()?, Kind::Match);
        assert_eq!("I".parse::<Kind>()?, Kind::Insertion);
        assert_eq!("D".parse::<Kind>()?, Kind::Deletion);
        assert_eq!("N".parse::<Kind>()?, Kind::Skip);
        assert_eq!("S".parse::<Kind>()?, Kind::SoftClip);
        assert_eq!("H".parse::<Kind>()?, Kind::HardClip);
        assert_eq!("P".parse::<Kind>()?, Kind::Pad);
        assert_eq!("=".parse::<Kind>()?, Kind::SeqMatch);
        assert_eq!("X".parse::<Kind>()?, Kind::SeqMismatch);

        assert!("".parse::<Kind>().is_err());
        assert!("O".parse::<Kind>().is_err());

        Ok(())
    }
}
