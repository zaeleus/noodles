bitflags::bitflags! {
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
    pub struct Flags: u8 {
        const RESERVED = 0x01;
        const DO_DEDUP = 0x02;
        const DO_LEN = 0x04;
        const DO_SEL = 0x08;
        const HAVE_QMAP = 0x10;
        const HAVE_PTAB = 0x20;
        const HAVE_DTAB = 0x40;
        const HAVE_QTAB = 0x80;
    }
}

impl Flags {
    pub fn has_duplicates(&self) -> bool {
        self.contains(Self::DO_DEDUP)
    }

    pub fn is_fixed_length(&self) -> bool {
        self.contains(Self::DO_LEN)
    }

    pub fn has_selector(&self) -> bool {
        self.contains(Self::DO_SEL)
    }

    pub fn has_quality_map(&self) -> bool {
        self.contains(Self::HAVE_QMAP)
    }

    pub fn has_positions_table(&self) -> bool {
        self.contains(Self::HAVE_PTAB)
    }

    pub fn has_deltas_table(&self) -> bool {
        self.contains(Self::HAVE_DTAB)
    }

    pub fn has_qualities_table(&self) -> bool {
        self.contains(Self::HAVE_QTAB)
    }
}

impl From<u8> for Flags {
    fn from(n: u8) -> Self {
        Self::from_bits_truncate(n)
    }
}

impl From<Flags> for u8 {
    fn from(flags: Flags) -> Self {
        flags.bits()
    }
}
