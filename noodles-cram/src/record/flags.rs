bitflags::bitflags! {
    /// CRAM record flags.
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
    pub struct Flags: u8 {
        /// The per-base quality scores are stored as an array, as opposed to read features
        /// (`0x01`).
        const QUALITY_SCORES_STORED_AS_ARRAY = 0x01;
        /// The mate is in another slice (`0x02`).
        const DETACHED = 0x02;
        /// The mate is after the current record (`0x04`).
        const HAS_MATE_DOWNSTREAM = 0x04;
        /// The sequence is unknown (`0x08`).
        const DECODE_SEQUENCE_AS_UNKNOWN = 0x08;
    }
}

impl Flags {
    /// Returns whether the `QUALITY_SCORES_STORED_AS_ARRAY` flag is set.
    pub fn are_quality_scores_stored_as_array(self) -> bool {
        self.contains(Self::QUALITY_SCORES_STORED_AS_ARRAY)
    }

    /// Returns whether the `DETACHED` flag is set.
    pub fn is_detached(self) -> bool {
        self.contains(Self::DETACHED)
    }

    /// Returns whether the `HAS_MATE_DOWNSTREAM` flag is set.
    pub fn has_mate_downstream(self) -> bool {
        self.contains(Self::HAS_MATE_DOWNSTREAM)
    }

    /// Returns whether the `DECODE_SEQUENCE_AS_UNKNOWN` flag is set.
    pub fn decode_sequence_as_unknown(self) -> bool {
        self.contains(Self::DECODE_SEQUENCE_AS_UNKNOWN)
    }
}

impl From<u8> for Flags {
    fn from(value: u8) -> Self {
        Self::from_bits_truncate(value)
    }
}

impl From<Flags> for u8 {
    fn from(flags: Flags) -> Self {
        flags.bits()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let flags = Flags::default();

        assert!(flags.is_empty());

        assert!(!flags.are_quality_scores_stored_as_array());
        assert!(!flags.is_detached());
        assert!(!flags.has_mate_downstream());
        assert!(!flags.decode_sequence_as_unknown());
    }

    #[test]
    fn test_contains() {
        assert!(Flags::QUALITY_SCORES_STORED_AS_ARRAY.are_quality_scores_stored_as_array());
        assert!(Flags::DETACHED.is_detached());
        assert!(Flags::HAS_MATE_DOWNSTREAM.has_mate_downstream());
        assert!(Flags::DECODE_SEQUENCE_AS_UNKNOWN.decode_sequence_as_unknown());
    }

    #[test]
    fn test_from_u8_for_flags() {
        assert_eq!(Flags::from(0x01), Flags::QUALITY_SCORES_STORED_AS_ARRAY);
        assert_eq!(Flags::from(0x02), Flags::DETACHED);
        assert_eq!(Flags::from(0x04), Flags::HAS_MATE_DOWNSTREAM);
        assert_eq!(Flags::from(0x08), Flags::DECODE_SEQUENCE_AS_UNKNOWN);
    }

    #[test]
    fn test_from_flags_for_u8() {
        assert_eq!(u8::from(Flags::from(0x01)), 0x01);
        assert_eq!(u8::from(Flags::from(0x02)), 0x02);
        assert_eq!(u8::from(Flags::from(0x04)), 0x04);
        assert_eq!(u8::from(Flags::from(0x08)), 0x08);
    }
}
