#[derive(Clone, Copy, Debug, Default)]
pub struct Flag(u16);

impl Flag {
    pub fn inner(self) -> u16 {
        self.0
    }

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

    #[deprecated(note = "Use Flag.is_duplicate instead.")]
    pub fn is_dup(self) -> bool {
        self.is_duplicate()
    }

    pub fn is_duplicate(self) -> bool {
        self.0 & 0x0400 != 0
    }

    pub fn is_supplementary(self) -> bool {
        self.0 & 0x0800 != 0
    }
}

impl From<u16> for Flag {
    fn from(value: u16) -> Self {
        Flag(value)
    }
}

#[cfg(test)]
mod tests {
    use super::Flag;

    #[test]
    fn test_empty_flag() {
        let flag = Flag::default();

        assert!(!flag.is_paired());
        assert!(!flag.is_proper_pair());
        assert!(!flag.is_unmapped());
        assert!(!flag.is_mate_unmapped());
        assert!(!flag.is_reverse());
        assert!(!flag.is_mate_reverse());
        assert!(!flag.is_read_1());
        assert!(!flag.is_read_2());
        assert!(!flag.is_secondary());
        assert!(!flag.is_qc_fail());
        assert!(!flag.is_duplicate());
        assert!(!flag.is_supplementary());
    }

    #[test]
    fn test_flags() {
        assert!(Flag::from(0x01).is_paired());
        assert!(Flag::from(0x02).is_proper_pair());
        assert!(Flag::from(0x04).is_unmapped());
        assert!(Flag::from(0x08).is_mate_unmapped());
        assert!(Flag::from(0x10).is_reverse());
        assert!(Flag::from(0x20).is_mate_reverse());
        assert!(Flag::from(0x40).is_read_1());
        assert!(Flag::from(0x80).is_read_2());
        assert!(Flag::from(0x0100).is_secondary());
        assert!(Flag::from(0x0200).is_qc_fail());
        assert!(Flag::from(0x0400).is_duplicate());
        assert!(Flag::from(0x0800).is_supplementary());
    }
}
