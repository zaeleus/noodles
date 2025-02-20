//! BGZF virtual position.

use std::{error, fmt};

pub(crate) const MAX_COMPRESSED_POSITION: u64 = (1 << 48) - 1;
pub(crate) const MAX_UNCOMPRESSED_POSITION: u16 = u16::MAX;

const COMPRESSED_POSITION_SHIFT: u64 = 16;
const UNCOMPRESSED_POSITION_MASK: u64 = 0xffff;

/// A BGZF virtual position.
///
/// A virtual position is a 64-bit unsigned integer representing both the position in the
/// compressed stream and position in the uncompressed block data. The compressed position is
/// typically at the start of a block.
///
/// The compressed position is the first six most significant bytes; and the uncompressed position,
/// the last two least significant bytes. For example, for the virtual position
/// 10253313912875616487:
///
/// ```text
///                       compressed position
///                        |               |
/// 10253313912875616487 = 8e 4b 16 ad eb 85 88 e7
///                                          |   |
///                                  uncompressed position
/// ```
///
/// The compressed position is at 156453154188165 (`8e 4b 16 ad eb 85`); and the uncompressed
/// position, 35047 (`88 e7`).
///
/// This is also called a virtual file offset; or, simply, a virtual offset.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct VirtualPosition(u64);

impl VirtualPosition {
    /// The minimum value of a virtual position.
    pub const MIN: Self = Self(u64::MIN);

    /// The maximum value of a virtual position.
    pub const MAX: Self = Self(u64::MAX);

    /// Creates a virtual position if the compressed position is valid.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::VirtualPosition;
    /// assert_eq!(VirtualPosition::new(0, 0), Some(VirtualPosition::MIN));
    /// assert!(VirtualPosition::new(1 << 48, 0).is_none());
    /// ```
    pub const fn new(compressed_pos: u64, uncompressed_pos: u16) -> Option<Self> {
        if compressed_pos <= MAX_COMPRESSED_POSITION {
            // SAFETY: 0 <= `uncompressed_pos` <= `u64`
            let virtual_pos =
                (compressed_pos << COMPRESSED_POSITION_SHIFT) | uncompressed_pos as u64;

            Some(Self(virtual_pos))
        } else {
            None
        }
    }

    /// The position in the compressed BGZF stream.
    ///
    /// This is typically at the start of a block.
    ///
    /// The maximum value of a compressed position is 281474976710655 (2^48 - 1).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let virtual_position = bgzf::VirtualPosition::from(3741638);
    /// assert_eq!(virtual_position.compressed(), 57);
    /// ```
    pub const fn compressed(self) -> u64 {
        self.0 >> COMPRESSED_POSITION_SHIFT
    }

    /// The position in the uncompressed block data.
    ///
    /// The maximum value of an uncompressed position is 65535 (2^16 - 1).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let virtual_position = bgzf::VirtualPosition::from(3741638);
    /// assert_eq!(virtual_position.uncompressed(), 6086);
    /// ```
    pub const fn uncompressed(self) -> u16 {
        (self.0 & UNCOMPRESSED_POSITION_MASK) as u16
    }
}

impl From<u64> for VirtualPosition {
    fn from(pos: u64) -> Self {
        Self(pos)
    }
}

/// An error returned when converting a (u64, u16) to a virtual position fails.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromU64U16TupleError {
    /// The compressed position is larger than 2^48 - 1.
    CompressedPositionOverflow,
}

impl error::Error for TryFromU64U16TupleError {}

impl fmt::Display for TryFromU64U16TupleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CompressedPositionOverflow => {
                f.write_str("the compressed position is larger than 2^48 - 1")
            }
        }
    }
}

impl TryFrom<(u64, u16)> for VirtualPosition {
    type Error = TryFromU64U16TupleError;

    /// Converts a `(compressed position, uncompressed position)` tuple to a virtual position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let virtual_position = bgzf::VirtualPosition::try_from((57, 6086));
    /// assert_eq!(virtual_position, Ok(bgzf::VirtualPosition::from(3741638)));
    /// ```
    fn try_from(pos: (u64, u16)) -> Result<Self, Self::Error> {
        let (compressed_pos, uncompressed_pos) = pos;

        if compressed_pos > MAX_COMPRESSED_POSITION {
            return Err(TryFromU64U16TupleError::CompressedPositionOverflow);
        }

        Ok(Self(
            (compressed_pos << COMPRESSED_POSITION_SHIFT) | u64::from(uncompressed_pos),
        ))
    }
}

impl From<VirtualPosition> for u64 {
    fn from(pos: VirtualPosition) -> Self {
        pos.0
    }
}

impl From<VirtualPosition> for (u64, u16) {
    fn from(pos: VirtualPosition) -> Self {
        (pos.compressed(), pos.uncompressed())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_u64_for_virtual_position() {
        let pos = VirtualPosition::from(88384945211);
        assert_eq!(pos.compressed(), 1348647);
        assert_eq!(pos.uncompressed(), 15419);

        let pos = VirtualPosition::from(188049630896);
        assert_eq!(pos.compressed(), 2869409);
        assert_eq!(pos.uncompressed(), 42672);

        let pos = VirtualPosition::from(26155658182977);
        assert_eq!(pos.compressed(), 399103671);
        assert_eq!(pos.uncompressed(), 321);
    }

    #[test]
    fn test_try_from_u64_u16_tuple_for_virtual_position() {
        assert_eq!(
            VirtualPosition::try_from((1348647, 15419)),
            Ok(VirtualPosition::from(88384945211))
        );

        assert_eq!(
            VirtualPosition::try_from((2869409, 42672)),
            Ok(VirtualPosition::from(188049630896))
        );

        assert_eq!(
            VirtualPosition::try_from((399103671, 321)),
            Ok(VirtualPosition::from(26155658182977))
        );

        assert_eq!(
            VirtualPosition::try_from((281474976710656, 0)),
            Err(TryFromU64U16TupleError::CompressedPositionOverflow)
        );
    }

    #[test]
    fn test_from_virtual_position_for_u64() {
        assert_eq!(u64::from(VirtualPosition::from(88384945211)), 88384945211);
        assert_eq!(u64::from(VirtualPosition::from(188049630896)), 188049630896);
        assert_eq!(
            u64::from(VirtualPosition::from(26155658182977)),
            26155658182977
        );
    }

    #[test]
    fn test_from_virtual_position_for_u64_u16_tuple() {
        assert_eq!(
            <(u64, u16)>::from(VirtualPosition::from(88384945211)),
            (1348647, 15419)
        );

        assert_eq!(
            <(u64, u16)>::from(VirtualPosition::from(188049630896)),
            (2869409, 42672)
        );

        assert_eq!(
            <(u64, u16)>::from(VirtualPosition::from(26155658182977)),
            (399103671, 321)
        );
    }
}
