bitflags::bitflags! {
    #[derive(Default)]
    pub struct Flags: u16 {
        const PAIRED = 0x01;
        const PROPER_PAIR = 0x02;
        const UNMAPPED = 0x04;
        const MATE_UNMAPPED = 0x08;
        const REVERSE = 0x10;
        const MATE_REVERSE = 0x20;
        const READ_1 = 0x40;
        const READ_2 = 0x80;
        const SECONDARY = 0x0100;
        const QC_FAIL = 0x0200;
        const DUPLICATE = 0x0400;
        const SUPPLEMENTARY = 0x0800;
    }
}

impl Flags {
    pub fn is_paired(self) -> bool {
        self.contains(Self::PAIRED)
    }

    pub fn is_proper_pair(self) -> bool {
        self.contains(Self::PROPER_PAIR)
    }

    pub fn is_unmapped(self) -> bool {
        self.contains(Self::UNMAPPED)
    }

    pub fn is_mate_unmapped(self) -> bool {
        self.contains(Self::MATE_UNMAPPED)
    }

    pub fn is_reverse(self) -> bool {
        self.contains(Self::REVERSE)
    }

    pub fn is_mate_reverse(self) -> bool {
        self.contains(Self::MATE_REVERSE)
    }

    pub fn is_read_1(self) -> bool {
        self.contains(Self::READ_1)
    }

    pub fn is_read_2(self) -> bool {
        self.contains(Self::READ_2)
    }

    pub fn is_secondary(self) -> bool {
        self.contains(Self::SECONDARY)
    }

    pub fn is_qc_fail(self) -> bool {
        self.contains(Self::QC_FAIL)
    }

    pub fn is_duplicate(self) -> bool {
        self.contains(Self::DUPLICATE)
    }

    pub fn is_supplementary(self) -> bool {
        self.contains(Self::SUPPLEMENTARY)
    }
}

impl From<u16> for Flags {
    fn from(value: u16) -> Self {
        Self::from_bits_truncate(value)
    }
}

impl From<Flags> for u16 {
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
    fn test_contains() {
        assert!(Flags::PAIRED.is_paired());
        assert!(Flags::PROPER_PAIR.is_proper_pair());
        assert!(Flags::UNMAPPED.is_unmapped());
        assert!(Flags::MATE_UNMAPPED.is_mate_unmapped());
        assert!(Flags::REVERSE.is_reverse());
        assert!(Flags::MATE_REVERSE.is_mate_reverse());
        assert!(Flags::READ_1.is_read_1());
        assert!(Flags::READ_2.is_read_2());
        assert!(Flags::SECONDARY.is_secondary());
        assert!(Flags::QC_FAIL.is_qc_fail());
        assert!(Flags::DUPLICATE.is_duplicate());
        assert!(Flags::SUPPLEMENTARY.is_supplementary());
    }

    #[test]
    fn test_from_u16_for_flags() {
        assert_eq!(Flags::from(0x40), Flags::READ_1);
    }

    #[test]
    fn test_from_flags_for_u16() {
        assert_eq!(u16::from(Flags::READ_1), 0x40);
    }
}
