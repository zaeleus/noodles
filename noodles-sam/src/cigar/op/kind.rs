use std::{fmt, str::FromStr};

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

impl FromStr for Kind {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "M" => Ok(Self::Match),
            "I" => Ok(Self::Insertion),
            "D" => Ok(Self::Deletion),
            "N" => Ok(Self::Skip),
            "S" => Ok(Self::SoftClip),
            "H" => Ok(Self::HardClip),
            "P" => Ok(Self::Pad),
            "=" => Ok(Self::SeqMismatch),
            "X" => Ok(Self::SeqMismatch),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Kind;

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
    }
}
