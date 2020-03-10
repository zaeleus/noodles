#[derive(Clone, Copy, Debug, Default)]
pub struct Flags(u16);

impl Flags {
    pub fn is_paired(self) -> bool {
        self.0 & 0x01 != 0
    }

    pub fn is_proper_pair(self) -> bool {
        self.0 & 0x02 != 0
    }

    pub fn is_unmapped(self) -> bool {
        self.0 & 0x04 != 0
    }

    pub fn is_mate_unmapped(self) -> bool {
        self.0 & 0x08 != 0
    }

    pub fn is_reverse(self) -> bool {
        self.0 & 0x10 != 0
    }

    pub fn is_mate_reverse(self) -> bool {
        self.0 & 0x20 != 0
    }

    pub fn is_read_1(self) -> bool {
        self.0 & 0x40 != 0
    }

    pub fn is_read_2(self) -> bool {
        self.0 & 0x80 != 0
    }

    pub fn is_secondary(self) -> bool {
        self.0 & 0x0100 != 0
    }

    pub fn is_qc_fail(self) -> bool {
        self.0 & 0x0200 != 0
    }

    pub fn is_duplicate(self) -> bool {
        self.0 & 0x0400 != 0
    }

    pub fn is_supplementary(self) -> bool {
        self.0 & 0x0800 != 0
    }
}

impl From<u16> for Flags {
    fn from(value: u16) -> Self {
        Self(value)
    }
}

impl From<Flags> for u16 {
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

        assert!(!flags.is_paired());
        assert!(!flags.is_proper_pair());
        assert!(!flags.is_unmapped());
        assert!(!flags.is_mate_unmapped());
        assert!(!flags.is_reverse());
        assert!(!flags.is_mate_reverse());
        assert!(!flags.is_read_1());
        assert!(!flags.is_read_2());
        assert!(!flags.is_secondary());
        assert!(!flags.is_qc_fail());
        assert!(!flags.is_duplicate());
        assert!(!flags.is_supplementary());
    }

    #[test]
    fn test_flags() {
        assert!(Flags::from(0x01).is_paired());
        assert!(Flags::from(0x02).is_proper_pair());
        assert!(Flags::from(0x04).is_unmapped());
        assert!(Flags::from(0x08).is_mate_unmapped());
        assert!(Flags::from(0x10).is_reverse());
        assert!(Flags::from(0x20).is_mate_reverse());
        assert!(Flags::from(0x40).is_read_1());
        assert!(Flags::from(0x80).is_read_2());
        assert!(Flags::from(0x0100).is_secondary());
        assert!(Flags::from(0x0200).is_qc_fail());
        assert!(Flags::from(0x0400).is_duplicate());
        assert!(Flags::from(0x0800).is_supplementary());
    }

    #[test]
    fn test_from_flags_for_u16() {
        assert_eq!(u16::from(Flags::from(0x40)), 0x40);
    }
}
