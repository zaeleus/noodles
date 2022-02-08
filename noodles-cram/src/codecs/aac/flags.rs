bitflags::bitflags! {
    pub struct Flags: u8 {
        const ORDER = 0x01;
        const RESERVED = 0x02;
        const EXT = 0x04;
        const STRIPE = 0x08;
        const NO_SIZE = 0x10;
        const CAT = 0x20;
        const RLE = 0x40;
        const PACK = 0x80;
    }
}

impl From<u8> for Flags {
    fn from(n: u8) -> Self {
        Self::from_bits_truncate(n)
    }
}
