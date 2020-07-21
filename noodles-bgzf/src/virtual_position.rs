pub(crate) const MAX_UNCOMPRESSED_POSITION: u16 = u16::MAX; // 2^16 - 1;

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
/// The compressed position is at 156453154188165 (`8e 4b 16 ad eb 85`), and uncompressed position
/// is at 35047 (`88 e7`).
///
/// This is also called a virtual file offset; or, simply, a virtual offset.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct VirtualPosition(u64);

impl VirtualPosition {
    /// The position in the compressed BGZF stream.
    ///
    /// This is typically at the start of a block.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let virtual_position = bgzf::VirtualPosition::from(3741638);
    /// assert_eq!(virtual_position.compressed(), 57);
    /// ```
    pub fn compressed(self) -> u64 {
        self.0 >> 16
    }

    /// The position in the uncompressed block data.
    ///
    /// The maximum value of an uncompressed position is 65535 (2^16-1).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let virtual_position = bgzf::VirtualPosition::from(3741638);
    /// assert_eq!(virtual_position.uncompressed(), 6086);
    /// ```
    pub fn uncompressed(self) -> u16 {
        (self.0 & 0xffff) as u16
    }
}

impl From<u64> for VirtualPosition {
    fn from(pos: u64) -> Self {
        Self(pos)
    }
}

impl From<(u64, u16)> for VirtualPosition {
    /// Converts a `(compressed position, uncompressed position)` tuple to a virtual position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let virtual_position = bgzf::VirtualPosition::from((57, 6086));
    /// assert_eq!(virtual_position, bgzf::VirtualPosition::from(3741638));
    /// ```
    fn from(pos: (u64, u16)) -> Self {
        let compressed_pos = pos.0 << 16;
        let uncompressed_pos = pos.1 as u64;
        Self(compressed_pos | uncompressed_pos)
    }
}

impl From<VirtualPosition> for u64 {
    fn from(pos: VirtualPosition) -> Self {
        pos.0
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
    fn test_from_u64_u16_tuple_for_virtual_position() {
        assert_eq!(
            VirtualPosition::from((1348647, 15419)),
            VirtualPosition::from(88384945211)
        );

        assert_eq!(
            VirtualPosition::from((2869409, 42672)),
            VirtualPosition::from(188049630896)
        );

        assert_eq!(
            VirtualPosition::from((399103671, 321)),
            VirtualPosition::from(26155658182977)
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
}
