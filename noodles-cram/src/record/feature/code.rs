//! CRAM record feature kind.

use std::{error, fmt};

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

/// An error returned when a byte fails to convert to a feature kind.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromByteError(u8);

impl error::Error for TryFromByteError {}

impl fmt::Display for TryFromByteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid code: {:#x}", self.0)
    }
}

/// An error returned when a character fails to convert to a feature kind.
#[deprecated(since = "0.13.0", note = "Use `TryFromByteError` instead.")]
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromCharError(char);

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

#[allow(deprecated)]
impl error::Error for TryFromCharError {}

#[allow(deprecated)]
impl fmt::Display for TryFromCharError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "expected {{b, q, B, X, I, D, i, Q, N, S, P, H}}, got {}",
            self.0
        )
    }
}

#[allow(deprecated)]
impl TryFrom<char> for Code {
    type Error = TryFromCharError;

    #[allow(useless_deprecated)]
    #[deprecated(since = "0.13.0", note = "Use `TryFrom<u8>` instead.")]
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

impl From<Code> for char {
    #[allow(useless_deprecated)]
    #[deprecated(since = "0.13.0", note = "Convert to a `u8` instead.")]
    fn from(code: Code) -> Self {
        match code {
            Code::Bases => 'b',
            Code::Scores => 'q',
            Code::ReadBase => 'B',
            Code::Substitution => 'X',
            Code::Insertion => 'I',
            Code::Deletion => 'D',
            Code::InsertBase => 'i',
            Code::QualityScore => 'Q',
            Code::ReferenceSkip => 'N',
            Code::SoftClip => 'S',
            Code::Padding => 'P',
            Code::HardClip => 'H',
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

    #[allow(deprecated)]
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

    #[allow(deprecated)]
    #[test]
    fn test_from_code_for_char() {
        assert_eq!(char::from(Code::Bases), 'b');
        assert_eq!(char::from(Code::Scores), 'q');
        assert_eq!(char::from(Code::ReadBase), 'B');
        assert_eq!(char::from(Code::Substitution), 'X');
        assert_eq!(char::from(Code::Insertion), 'I');
        assert_eq!(char::from(Code::Deletion), 'D');
        assert_eq!(char::from(Code::InsertBase), 'i');
        assert_eq!(char::from(Code::QualityScore), 'Q');
        assert_eq!(char::from(Code::ReferenceSkip), 'N');
        assert_eq!(char::from(Code::SoftClip), 'S');
        assert_eq!(char::from(Code::Padding), 'P');
        assert_eq!(char::from(Code::HardClip), 'H');
    }
}
