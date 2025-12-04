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

impl Flags {
    pub(super) fn order(&self) -> u8 {
        if self.contains(Self::ORDER) { 1 } else { 0 }
    }

    pub(super) fn uses_external_codec(&self) -> bool {
        self.contains(Self::EXT)
    }

    pub(super) fn is_striped(&self) -> bool {
        self.contains(Self::STRIPE)
    }

    pub(super) fn has_uncompressed_size(&self) -> bool {
        !self.contains(Self::NO_SIZE)
    }

    pub(super) fn is_uncompressed(&self) -> bool {
        self.contains(Self::CAT)
    }

    pub(super) fn is_rle(&self) -> bool {
        self.contains(Self::RLE)
    }

    pub(super) fn is_bit_packed(&self) -> bool {
        self.contains(Self::PACK)
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
