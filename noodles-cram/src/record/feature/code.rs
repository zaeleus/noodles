use std::{convert::TryFrom, error, fmt};

/// A CRAM record feature kind.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Code {
    /// A stretch of bases (`b`).
    Bases,
    /// A stretch of quality scores (`q`).
    Scores,
    /// A (base, quality score) pair (`B`).
    ReadBase,
    /// A base substitution (CIGAR op `X`, `M` and `=`; `X`).
    Substitution,
    /// Inserted bases (CIGAR op `I`; `I`).
    Insertion,
    /// A number of deleted bases (CIGAR op `D`; `D`).
    Deletion,
    /// A single inserted base (CIGAR op `I`; `i`).
    InsertBase,
    /// A single quality score (`Q`).
    QualityScore,
    /// A number of skipped bases (CIGAR op `N`; `N`).
    ReferenceSkip,
    /// Soft clipped bases (CIGAR op `S`; `S`).
    SoftClip,
    /// A number of padded bases (CIGAR op `P`; `P`).
    Padding,
    /// A number of hard clipped bases (CIGAR op `H`, `H`).
    HardClip,
}

/// An error returned when a character fails to convert to a feature kind.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromCharError(char);

impl error::Error for TryFromCharError {}

impl fmt::Display for TryFromCharError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "expected {{b, q, B, X, I, D, i, Q, N, S, P, H}}, got {}",
            self.0
        )
    }
}

impl TryFrom<char> for Code {
    type Error = TryFromCharError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            'b' => Ok(Self::Bases),
            'q' => Ok(Self::Scores),
            'B' => Ok(Self::ReadBase),
            'X' => Ok(Self::Substitution),
            'I' => Ok(Self::Insertion),
            'D' => Ok(Self::Deletion),
            'i' => Ok(Self::InsertBase),
            'Q' => Ok(Self::QualityScore),
            'N' => Ok(Self::ReferenceSkip),
            'S' => Ok(Self::SoftClip),
            'P' => Ok(Self::Padding),
            'H' => Ok(Self::HardClip),
            _ => Err(TryFromCharError(c)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_char_for_code() {
        assert_eq!(Code::try_from('b'), Ok(Code::Bases));
        assert_eq!(Code::try_from('q'), Ok(Code::Scores));
        assert_eq!(Code::try_from('B'), Ok(Code::ReadBase));
        assert_eq!(Code::try_from('X'), Ok(Code::Substitution));
        assert_eq!(Code::try_from('I'), Ok(Code::Insertion));
        assert_eq!(Code::try_from('D'), Ok(Code::Deletion));
        assert_eq!(Code::try_from('i'), Ok(Code::InsertBase));
        assert_eq!(Code::try_from('Q'), Ok(Code::QualityScore));
        assert_eq!(Code::try_from('N'), Ok(Code::ReferenceSkip));
        assert_eq!(Code::try_from('S'), Ok(Code::SoftClip));
        assert_eq!(Code::try_from('P'), Ok(Code::Padding));
        assert_eq!(Code::try_from('H'), Ok(Code::HardClip));
        assert_eq!(Code::try_from('Z'), Err(TryFromCharError('Z')));
    }
}
