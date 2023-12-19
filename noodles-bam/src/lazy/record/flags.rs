use std::io;

use noodles_sam as sam;

/// Raw BAM record flags.
#[derive(Debug, Eq, PartialEq)]
pub struct Flags(u16);

impl Flags {
    pub(super) fn new(n: u16) -> Self {
        Self(n)
    }
}

impl sam::alignment::record::Flags for Flags {
    fn try_to_u16(&self) -> io::Result<u16> {
        Ok(self.0)
    }
}

impl From<Flags> for u16 {
    fn from(flags: Flags) -> Self {
        flags.0
    }
}

impl From<Flags> for sam::record::Flags {
    fn from(flags: Flags) -> Self {
        Self::from(flags.0)
    }
}
