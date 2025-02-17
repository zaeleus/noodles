bitflags::bitflags! {
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
    pub struct Flags: u8 {
        const QUALITY_SCORES_ARE_STORED_AS_ARRAY = 0x01;
        const IS_DETACHED = 0x02;
        const MATE_IS_DOWNSTREAM = 0x04;
        const SEQUENCE_IS_MISSING = 0x08;
    }
}

impl Flags {
    pub fn quality_scores_are_stored_as_array(self) -> bool {
        self.contains(Self::QUALITY_SCORES_ARE_STORED_AS_ARRAY)
    }

    pub fn is_detached(self) -> bool {
        self.contains(Self::IS_DETACHED)
    }

    pub fn mate_is_downstream(self) -> bool {
        self.contains(Self::MATE_IS_DOWNSTREAM)
    }

    pub fn sequence_is_missing(self) -> bool {
        self.contains(Self::SEQUENCE_IS_MISSING)
    }
}

impl From<u8> for Flags {
    fn from(value: u8) -> Self {
        Self::from_bits_truncate(value)
    }
}

impl From<Flags> for u8 {
    fn from(flags: Flags) -> Self {
        flags.bits()
    }
}
