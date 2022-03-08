//! BAM record CIGAR operation.

use std::{error, fmt, mem};

use byteorder::{ByteOrder, LittleEndian};
use noodles_sam::{self as sam, record::cigar::op::Kind};

const MAX_LENGTH: u32 = (1 << 28) - 1;

/// A BAM record CIGAR operation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Op {
    kind: Kind,
    len: u32,
}

/// An error returned when the operation length is > 2^28.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LengthError(u32);

impl error::Error for LengthError {}

impl fmt::Display for LengthError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "expected operation length to be <= {}, got {}",
            MAX_LENGTH, self.0
        )
    }
}

impl Op {
    /// Creates a CIGAR operation.
    ///
    /// # Errors
    ///
    /// [`LengthError`] is returned if the length is > 2^28 - 1.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::cigar::Op;
    /// use noodles_sam::record::cigar::op::Kind;
    /// let op = Op::new(Kind::Match, 13)?;
    /// # Ok::<(), noodles_bam::record::cigar::op::LengthError>(())
    /// ```
    pub fn new(kind: Kind, len: u32) -> Result<Self, LengthError> {
        if len <= MAX_LENGTH {
            Ok(Self { kind, len })
        } else {
            Err(LengthError(len))
        }
    }

    /// Returns the kind of the operation.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::cigar::Op;
    /// use noodles_sam::record::cigar::op::Kind;
    /// let op = Op::new(Kind::Match, 13)?;
    /// assert_eq!(op.kind(), Kind::Match);
    /// # Ok::<(), noodles_bam::record::cigar::op::LengthError>(())
    /// ```
    pub fn kind(self) -> Kind {
        self.kind
    }

    /// Returns the length of the operation.
    ///
    /// The length is guaranteed to be <= 2^28 - 1.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::cigar::Op;
    /// use noodles_sam::record::cigar::op::Kind;
    /// let op = Op::new(Kind::Match, 13)?;
    /// assert_eq!(op.len(), 13);
    /// # Ok::<(), noodles_bam::record::cigar::op::LengthError>(())
    /// ```
    pub fn len(self) -> u32 {
        self.len
    }

    /// Returns whether the operation is a no-op.
    ///
    /// That is, whether the operation has a length of 0.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::cigar::Op;
    /// use noodles_sam::record::cigar::op::Kind;
    ///
    /// let op = Op::new(Kind::Match, 0)?;
    /// assert!(op.is_empty());
    ///
    /// let op = Op::new(Kind::Match, 13)?;
    /// assert!(!op.is_empty());
    /// # Ok::<(), noodles_bam::record::cigar::op::LengthError>(())
    /// ```
    pub fn is_empty(self) -> bool {
        self.len == 0
    }
}

/// An error returned when a raw u32 fails to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromUIntError {
    /// The operation is invalid.
    InvalidOp(u32),
    /// The length is invalid.
    InvalidLength(LengthError),
}

impl error::Error for TryFromUIntError {}

impl fmt::Display for TryFromUIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidOp(n) => write!(f, "invalid op: expected 0..=8, got {}", n),
            Self::InvalidLength(e) => write!(f, "invalid length: {}", e),
        }
    }
}

impl TryFrom<u32> for Op {
    type Error = TryFromUIntError;

    fn try_from(u: u32) -> Result<Self, Self::Error> {
        let len = u >> 4;

        let kind = match u & 0x0f {
            0 => Kind::Match,
            1 => Kind::Insertion,
            2 => Kind::Deletion,
            3 => Kind::Skip,
            4 => Kind::SoftClip,
            5 => Kind::HardClip,
            6 => Kind::Pad,
            7 => Kind::SequenceMatch,
            8 => Kind::SequenceMismatch,
            n => return Err(TryFromUIntError::InvalidOp(n)),
        };

        Self::new(kind, len).map_err(TryFromUIntError::InvalidLength)
    }
}

/// An error returned when a raw byte slice fails to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromByteSliceError {
    /// The input is invalid.
    Invalid(usize),
    /// The value is invalid.
    InvalidUInt(TryFromUIntError),
}

impl error::Error for TryFromByteSliceError {}

impl fmt::Display for TryFromByteSliceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(len) => write!(f, "expected length to >= 4, got {}", len),
            Self::InvalidUInt(e) => write!(f, "invalid u32: {}", e),
        }
    }
}

impl TryFrom<&[u8]> for Op {
    type Error = TryFromByteSliceError;

    #[allow(useless_deprecated)]
    #[deprecated(since = "0.8.0", note = "Use `TryFrom<u32>` instead.")]
    fn try_from(buf: &[u8]) -> Result<Self, Self::Error> {
        if buf.len() < mem::size_of::<u32>() {
            return Err(TryFromByteSliceError::Invalid(buf.len()));
        }

        let u = LittleEndian::read_u32(buf);
        Self::try_from(u).map_err(TryFromByteSliceError::InvalidUInt)
    }
}

impl fmt::Display for Op {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.len(), self.kind())
    }
}

impl From<Op> for u32 {
    fn from(op: Op) -> Self {
        let i = op.kind() as u32;
        op.len() << 4 | i
    }
}

impl From<Op> for sam::record::cigar::Op {
    fn from(op: Op) -> Self {
        Self::new(op.kind(), op.len())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_u32() -> Result<(), LengthError> {
        assert_eq!(Op::try_from(1 << 4), Ok(Op::new(Kind::Match, 1)?));
        assert_eq!(Op::try_from(2 << 4 | 1), Ok(Op::new(Kind::Insertion, 2)?));
        assert_eq!(Op::try_from(3 << 4 | 2), Ok(Op::new(Kind::Deletion, 3)?));
        assert_eq!(Op::try_from(4 << 4 | 3), Ok(Op::new(Kind::Skip, 4)?));
        assert_eq!(Op::try_from(5 << 4 | 4), Ok(Op::new(Kind::SoftClip, 5)?));
        assert_eq!(Op::try_from(6 << 4 | 5), Ok(Op::new(Kind::HardClip, 6)?));
        assert_eq!(Op::try_from(7 << 4 | 6), Ok(Op::new(Kind::Pad, 7)?));
        assert_eq!(
            Op::try_from(8 << 4 | 7),
            Ok(Op::new(Kind::SequenceMatch, 8)?)
        );
        assert_eq!(
            Op::try_from(9 << 4 | 8),
            Ok(Op::new(Kind::SequenceMismatch, 9)?)
        );

        assert_eq!(
            Op::try_from(10 << 4 | 9),
            Err(TryFromUIntError::InvalidOp(9))
        );

        Ok(())
    }

    #[test]
    fn test_try_from_u8_slice() -> Result<(), LengthError> {
        let buf = [0x40, 0x02, 0x00, 0x00];
        assert_eq!(Op::try_from(&buf[..]), Ok(Op::new(Kind::Match, 36)?));

        let buf = [0x40, 0x02, 0x00, 0x00, 0x84, 0x00, 0x00, 0x00];
        assert_eq!(Op::try_from(&buf[..]), Ok(Op::new(Kind::Match, 36)?));

        let buf = [];
        assert_eq!(
            Op::try_from(&buf[..]),
            Err(TryFromByteSliceError::Invalid(0))
        );

        let buf = [0x49, 0x02, 0x00, 0x00];
        assert_eq!(
            Op::try_from(&buf[..]),
            Err(TryFromByteSliceError::InvalidUInt(
                TryFromUIntError::InvalidOp(9)
            ))
        );

        Ok(())
    }

    #[test]
    fn test_fmt() -> Result<(), LengthError> {
        assert_eq!(Op::new(Kind::Match, 32)?.to_string(), "32M");
        assert_eq!(Op::new(Kind::Deletion, 7)?.to_string(), "7D");
        assert_eq!(Op::new(Kind::Skip, 11)?.to_string(), "11N");
        assert_eq!(Op::new(Kind::Pad, 188)?.to_string(), "188P");
        Ok(())
    }

    #[test]
    fn test_from_op_for_u32() -> Result<(), LengthError> {
        assert_eq!(u32::from(Op::new(Kind::Match, 1)?), 1 << 4);
        assert_eq!(u32::from(Op::new(Kind::Insertion, 2)?), 2 << 4 | 1);
        assert_eq!(u32::from(Op::new(Kind::Deletion, 3)?), 3 << 4 | 2);
        assert_eq!(u32::from(Op::new(Kind::Skip, 4)?), 4 << 4 | 3);
        assert_eq!(u32::from(Op::new(Kind::SoftClip, 5)?), 5 << 4 | 4);
        assert_eq!(u32::from(Op::new(Kind::HardClip, 6)?), 6 << 4 | 5);
        assert_eq!(u32::from(Op::new(Kind::Pad, 7)?), 7 << 4 | 6);
        assert_eq!(u32::from(Op::new(Kind::SequenceMatch, 8)?), 8 << 4 | 7);
        assert_eq!(u32::from(Op::new(Kind::SequenceMismatch, 9)?), 9 << 4 | 8);
        Ok(())
    }

    #[test]
    fn test_from_op_for_sam_record_cigar_op() -> Result<(), LengthError> {
        let op = Op::new(Kind::Match, 36)?;

        let actual = sam::record::cigar::Op::from(op);
        let expected = sam::record::cigar::Op::new(op.kind(), op.len());

        assert_eq!(actual, expected);

        Ok(())
    }
}
