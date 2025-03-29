use std::{error, fmt};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Code {
    Bases,
    Scores,
    ReadBase,
    Substitution,
    Insertion,
    Deletion,
    InsertBase,
    QualityScore,
    ReferenceSkip,
    SoftClip,
    Padding,
    HardClip,
}

/// An error returned when a byte fails to convert to a feature kind.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromByteError(u8);

impl error::Error for TryFromByteError {}

impl fmt::Display for TryFromByteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid code: {:#x}", self.0)
    }
}

impl TryFrom<u8> for Code {
    type Error = TryFromByteError;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        match n {
            b'b' => Ok(Self::Bases),
            b'q' => Ok(Self::Scores),
            b'B' => Ok(Self::ReadBase),
            b'X' => Ok(Self::Substitution),
            b'I' => Ok(Self::Insertion),
            b'D' => Ok(Self::Deletion),
            b'i' => Ok(Self::InsertBase),
            b'Q' => Ok(Self::QualityScore),
            b'N' => Ok(Self::ReferenceSkip),
            b'S' => Ok(Self::SoftClip),
            b'P' => Ok(Self::Padding),
            b'H' => Ok(Self::HardClip),
            _ => Err(TryFromByteError(n)),
        }
    }
}

impl From<Code> for u8 {
    fn from(code: Code) -> Self {
        match code {
            Code::Bases => b'b',
            Code::Scores => b'q',
            Code::ReadBase => b'B',
            Code::Substitution => b'X',
            Code::Insertion => b'I',
            Code::Deletion => b'D',
            Code::InsertBase => b'i',
            Code::QualityScore => b'Q',
            Code::ReferenceSkip => b'N',
            Code::SoftClip => b'S',
            Code::Padding => b'P',
            Code::HardClip => b'H',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_u8_for_code() {
        assert_eq!(Code::try_from(b'b'), Ok(Code::Bases));
        assert_eq!(Code::try_from(b'q'), Ok(Code::Scores));
        assert_eq!(Code::try_from(b'B'), Ok(Code::ReadBase));
        assert_eq!(Code::try_from(b'X'), Ok(Code::Substitution));
        assert_eq!(Code::try_from(b'I'), Ok(Code::Insertion));
        assert_eq!(Code::try_from(b'D'), Ok(Code::Deletion));
        assert_eq!(Code::try_from(b'i'), Ok(Code::InsertBase));
        assert_eq!(Code::try_from(b'Q'), Ok(Code::QualityScore));
        assert_eq!(Code::try_from(b'N'), Ok(Code::ReferenceSkip));
        assert_eq!(Code::try_from(b'S'), Ok(Code::SoftClip));
        assert_eq!(Code::try_from(b'P'), Ok(Code::Padding));
        assert_eq!(Code::try_from(b'H'), Ok(Code::HardClip));
        assert_eq!(Code::try_from(b'Z'), Err(TryFromByteError(b'Z')));
    }

    #[test]
    fn test_from_code_for_u8() {
        assert_eq!(u8::from(Code::Bases), b'b');
        assert_eq!(u8::from(Code::Scores), b'q');
        assert_eq!(u8::from(Code::ReadBase), b'B');
        assert_eq!(u8::from(Code::Substitution), b'X');
        assert_eq!(u8::from(Code::Insertion), b'I');
        assert_eq!(u8::from(Code::Deletion), b'D');
        assert_eq!(u8::from(Code::InsertBase), b'i');
        assert_eq!(u8::from(Code::QualityScore), b'Q');
        assert_eq!(u8::from(Code::ReferenceSkip), b'N');
        assert_eq!(u8::from(Code::SoftClip), b'S');
        assert_eq!(u8::from(Code::Padding), b'P');
        assert_eq!(u8::from(Code::HardClip), b'H');
    }
}
