bitflags::bitflags! {
    /// CRAM record next mate flags.
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
    pub struct NextMateFlags: u8 {
        /// The mate is on the negative strand (`0x01`).
        const ON_NEGATIVE_STRAND = 0x01;
        /// The mate is unmapped (`0x02`).
        const UNMAPPED = 0x02;
    }
}

impl NextMateFlags {
    /// Returns whether the `ON_NEGATIVE_STRAND` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::record::NextMateFlags;
    /// assert!(NextMateFlags::ON_NEGATIVE_STRAND.is_on_negative_strand());
    /// ```
    pub fn is_on_negative_strand(self) -> bool {
        self.contains(Self::ON_NEGATIVE_STRAND)
    }

    /// Returns whether the `IS_UNMAPPED` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::record::NextMateFlags;
    /// assert!(NextMateFlags::UNMAPPED.is_unmapped());
    /// ```
    pub fn is_unmapped(self) -> bool {
        self.contains(Self::UNMAPPED)
    }
}

impl From<u8> for NextMateFlags {
    fn from(value: u8) -> Self {
        Self::from_bits_truncate(value)
    }
}

impl From<NextMateFlags> for u8 {
    fn from(flags: NextMateFlags) -> Self {
        flags.bits()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let flags = NextMateFlags::default();

        assert!(flags.is_empty());

        assert!(!flags.is_on_negative_strand());
        assert!(!flags.is_unmapped());
    }

    #[test]
    fn test_contains() {
        assert!(NextMateFlags::ON_NEGATIVE_STRAND.is_on_negative_strand());
        assert!(NextMateFlags::UNMAPPED.is_unmapped());
    }

    #[test]
    fn test_from_u8_for_flags() {
        assert_eq!(NextMateFlags::from(0x01), NextMateFlags::ON_NEGATIVE_STRAND);
        assert_eq!(NextMateFlags::from(0x02), NextMateFlags::UNMAPPED);
    }

    #[test]
    fn test_from_flags_for_u8() {
        assert_eq!(u8::from(NextMateFlags::from(0x01)), 0x01);
        assert_eq!(u8::from(NextMateFlags::from(0x02)), 0x02);
    }
}
