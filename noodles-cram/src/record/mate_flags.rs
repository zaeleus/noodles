bitflags::bitflags! {
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
    pub struct MateFlags: u8 {
        const IS_ON_NEGATIVE_STRAND = 0x01;
        const IS_UNMAPPED = 0x02;
    }
}

impl MateFlags {
    pub fn is_on_negative_strand(self) -> bool {
        self.contains(Self::IS_ON_NEGATIVE_STRAND)
    }

    pub fn is_unmapped(self) -> bool {
        self.contains(Self::IS_UNMAPPED)
    }
}

impl From<u8> for MateFlags {
    fn from(value: u8) -> Self {
        Self::from_bits_truncate(value)
    }
}

impl From<MateFlags> for u8 {
    fn from(flags: MateFlags) -> Self {
        flags.bits()
    }
}
