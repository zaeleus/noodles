#[derive(Clone, Copy, Debug, Default)]
pub struct Flags(u8);

impl Flags {
    pub fn are_quality_scores_stored_as_array(self) -> bool {
        self.0 & 0x01 != 0
    }

    pub fn is_detached(self) -> bool {
        self.0 & 0x02 != 0
    }

    pub fn has_mate_downstream(self) -> bool {
        self.0 & 0x04 != 0
    }

    pub fn decode_sequence_as_unknown(self) -> bool {
        self.0 & 0x08 != 0
    }
}

impl From<u8> for Flags {
    fn from(value: u8) -> Self {
        Self(value)
    }
}

impl From<Flags> for u8 {
    fn from(flags: Flags) -> Self {
        flags.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_flags() {
        let flags = Flags::default();

        assert!(!flags.are_quality_scores_stored_as_array());
        assert!(!flags.is_detached());
        assert!(!flags.has_mate_downstream());
        assert!(!flags.decode_sequence_as_unknown());
    }

    #[test]
    fn test_flags() {
        assert!(Flags::from(0x01).are_quality_scores_stored_as_array());
        assert!(Flags::from(0x02).is_detached());
        assert!(Flags::from(0x04).has_mate_downstream());
        assert!(Flags::from(0x08).decode_sequence_as_unknown());
    }

    #[test]
    fn test_from_flags_for_u16() {
        assert_eq!(u8::from(Flags::from(0x01)), 0x01);
        assert_eq!(u8::from(Flags::from(0x02)), 0x02);
        assert_eq!(u8::from(Flags::from(0x04)), 0x04);
        assert_eq!(u8::from(Flags::from(0x08)), 0x08);
    }
}
