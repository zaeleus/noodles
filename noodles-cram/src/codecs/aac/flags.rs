bitflags::bitflags! {
    /// Adaptive arithmetic coder flags.
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
    pub struct Flags: u8 {
        /// Order-1 entropy coding.
        const ORDER = 0x01;
        /// Reserved.
        const RESERVED = 0x02;
        /// Use an external codec (bzip2).
        const EXT = 0x04;
        /// Interleave 32 rANS states.
        const STRIPE = 0x08;
        /// Discard the uncompressed data size.
        const NO_SIZE = 0x10;
        /// Data is uncompressed.
        const CAT = 0x20;
        /// Use run-length encoding.
        const RLE = 0x40;
        /// Use bit packing.
        const PACK = 0x80;
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
