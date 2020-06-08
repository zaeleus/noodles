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

impl fmt::Display for Kind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", char::from(*self))
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
