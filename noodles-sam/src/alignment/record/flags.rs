bitflags::bitflags! {
    /// Alignment record flags.
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
    pub struct Flags: u16 {
        /// Read is segmented (`0x01`).
        const SEGMENTED = 0x01;
        /// Each segment in the read is properly aligned (`0x02`).
        const PROPERLY_SEGMENTED = 0x02;
        /// Read is unmapped (`Ox04`).
        const UNMAPPED = 0x04;
        /// The mate is unmapped (`0x08`).
        const MATE_UNMAPPED = 0x08;
        /// The sequence is reverse complemented (`0x10`).
        const REVERSE_COMPLEMENTED = 0x10;
        /// The sequence of the mate is reverse complemented (`0x20`).
        const MATE_REVERSE_COMPLEMENTED = 0x20;
        /// First segment in the read (`0x40`).
        const FIRST_SEGMENT = 0x40;
        /// Last segment in the read (`0x80`).
        const LAST_SEGMENT = 0x80;
        /// Secondary read (`0x0100`).
        const SECONDARY = 0x0100;
        /// Read failed quality checks (`0x0200`).
        const QC_FAIL = 0x0200;
        /// PCR or optical duplicate read (`0x0400`).
        const DUPLICATE = 0x0400;
        /// Supplementary alignment (`0x0800`).
        const SUPPLEMENTARY = 0x0800;
    }
}

impl Flags {
    /// Returns whether the `SEGMENTED` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::Flags;
    /// assert!(Flags::SEGMENTED.is_segmented());
    /// assert!(!Flags::UNMAPPED.is_segmented());
    /// ```
    pub fn is_segmented(self) -> bool {
        self.contains(Self::SEGMENTED)
    }

    /// Returns whether the `PROPERLY_SEGMENTED` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::Flags;
    /// assert!(Flags::PROPERLY_SEGMENTED.is_properly_segmented());
    /// assert!(!Flags::UNMAPPED.is_properly_segmented());
    /// ```
    pub fn is_properly_segmented(self) -> bool {
        self.contains(Self::PROPERLY_SEGMENTED)
    }

    /// Returns whether the `UNMAPPED` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::Flags;
    /// assert!(Flags::UNMAPPED.is_unmapped());
    /// assert!(!Flags::SEGMENTED.is_unmapped());
    /// ```
    pub fn is_unmapped(self) -> bool {
        self.contains(Self::UNMAPPED)
    }

    /// Returns whether the `MATE_UNMAPPED` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::Flags;
    /// assert!(Flags::MATE_UNMAPPED.is_mate_unmapped());
    /// assert!(!Flags::UNMAPPED.is_mate_unmapped());
    /// ```
    pub fn is_mate_unmapped(self) -> bool {
        self.contains(Self::MATE_UNMAPPED)
    }

    /// Returns whether the `REVERSE_COMPLEMENTED` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::Flags;
    /// assert!(Flags::REVERSE_COMPLEMENTED.is_reverse_complemented());
    /// assert!(!Flags::UNMAPPED.is_reverse_complemented());
    /// ```
    pub fn is_reverse_complemented(self) -> bool {
        self.contains(Self::REVERSE_COMPLEMENTED)
    }

    /// Returns whether the `MATE_REVERSE_COMPLEMENTED` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::Flags;
    /// assert!(Flags::MATE_REVERSE_COMPLEMENTED.is_mate_reverse_complemented());
    /// assert!(!Flags::UNMAPPED.is_mate_reverse_complemented());
    /// ```
    pub fn is_mate_reverse_complemented(self) -> bool {
        self.contains(Self::MATE_REVERSE_COMPLEMENTED)
    }

    /// Returns whether the `FIRST_SEGMENT` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::Flags;
    /// assert!(Flags::FIRST_SEGMENT.is_first_segment());
    /// assert!(!Flags::UNMAPPED.is_first_segment());
    /// ```
    pub fn is_first_segment(self) -> bool {
        self.contains(Self::FIRST_SEGMENT)
    }

    /// Returns whether the `LAST_SEGMENT` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::Flags;
    /// assert!(Flags::LAST_SEGMENT.is_last_segment());
    /// assert!(!Flags::UNMAPPED.is_last_segment());
    /// ```
    pub fn is_last_segment(self) -> bool {
        self.contains(Self::LAST_SEGMENT)
    }

    /// Returns whether the `SECONDARY` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::Flags;
    /// assert!(Flags::SECONDARY.is_secondary());
    /// assert!(!Flags::UNMAPPED.is_secondary());
    /// ```
    pub fn is_secondary(self) -> bool {
        self.contains(Self::SECONDARY)
    }

    /// Returns whether the `QC_FAIL` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::Flags;
    /// assert!(Flags::QC_FAIL.is_qc_fail());
    /// assert!(!Flags::UNMAPPED.is_qc_fail());
    /// ```
    pub fn is_qc_fail(self) -> bool {
        self.contains(Self::QC_FAIL)
    }

    /// Returns whether the `DUPLICATE` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::Flags;
    /// assert!(Flags::DUPLICATE.is_duplicate());
    /// assert!(!Flags::UNMAPPED.is_duplicate());
    /// ```
    pub fn is_duplicate(self) -> bool {
        self.contains(Self::DUPLICATE)
    }

    /// Returns whether the `SUPPLEMENTARY` flag is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::Flags;
    /// assert!(Flags::SUPPLEMENTARY.is_supplementary());
    /// assert!(!Flags::UNMAPPED.is_supplementary());
    /// ```
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
        assert!(!flags.is_segmented());
        assert!(!flags.is_properly_segmented());
        assert!(!flags.is_unmapped());
        assert!(!flags.is_mate_unmapped());
        assert!(!flags.is_reverse_complemented());
        assert!(!flags.is_mate_reverse_complemented());
        assert!(!flags.is_first_segment());
        assert!(!flags.is_last_segment());
        assert!(!flags.is_secondary());
        assert!(!flags.is_qc_fail());
        assert!(!flags.is_duplicate());
        assert!(!flags.is_supplementary());
    }

    #[test]
    fn test_contains() {
        assert!(Flags::SEGMENTED.is_segmented());
        assert!(Flags::PROPERLY_SEGMENTED.is_properly_segmented());
        assert!(Flags::UNMAPPED.is_unmapped());
        assert!(Flags::MATE_UNMAPPED.is_mate_unmapped());
        assert!(Flags::REVERSE_COMPLEMENTED.is_reverse_complemented());
        assert!(Flags::MATE_REVERSE_COMPLEMENTED.is_mate_reverse_complemented());
        assert!(Flags::FIRST_SEGMENT.is_first_segment());
        assert!(Flags::LAST_SEGMENT.is_last_segment());
        assert!(Flags::SECONDARY.is_secondary());
        assert!(Flags::QC_FAIL.is_qc_fail());
        assert!(Flags::DUPLICATE.is_duplicate());
        assert!(Flags::SUPPLEMENTARY.is_supplementary());
    }

    #[test]
    fn test_from_u16_for_flags() {
        assert_eq!(Flags::from(0x04), Flags::UNMAPPED);
    }

    #[test]
    fn test_from_flags_for_u16() {
        assert_eq!(u16::from(Flags::UNMAPPED), 0x04);
    }
}
