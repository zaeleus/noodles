use std::io;

use bitflags::bitflags;

use crate::io::reader::num::read_u8;

bitflags! {
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
    pub struct Flags: u8 {
        const MULTI_PARAM = 0x01;
        const HAVE_S_TAB = 0x02;
        const DO_REV = 0x04;
    }
}

impl Flags {
    pub(super) fn has_parameter_count(&self) -> bool {
        self.contains(Self::MULTI_PARAM)
    }

    pub(super) fn has_selector_table(&self) -> bool {
        self.contains(Self::HAVE_S_TAB)
    }

    pub fn has_reversed_values(&self) -> bool {
        self.contains(Self::DO_REV)
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

pub(super) fn read_flags(src: &mut &[u8]) -> io::Result<Flags> {
    read_u8(src).map(Flags::from)
}
