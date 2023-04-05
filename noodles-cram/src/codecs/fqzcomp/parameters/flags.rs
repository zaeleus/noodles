bitflags::bitflags! {
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
    pub struct Flags: u8 {
        const MULTI_PARAM = 0x01;
        const HAVE_S_TAB = 0x02;
        const DO_REV = 0x04;
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
