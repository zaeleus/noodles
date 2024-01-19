use std::io;

use crate::alignment::record::field::MappingQuality as _;

/// Raw SAM record mapping quality.
#[derive(Debug, Eq, PartialEq)]
pub struct MappingQuality<'a>(&'a [u8]);

impl<'a> MappingQuality<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }
}

impl<'a> crate::alignment::record::field::MappingQuality for MappingQuality<'a> {
    fn try_to_u8(&self) -> io::Result<u8> {
        lexical_core::parse(self.as_ref())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}

impl<'a> AsRef<[u8]> for MappingQuality<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryFrom<MappingQuality<'a>> for u8 {
    type Error = io::Error;

    fn try_from(raw_mapping_quality: MappingQuality<'a>) -> Result<Self, Self::Error> {
        raw_mapping_quality.try_to_u8()
    }
}

impl<'a> TryFrom<MappingQuality<'a>> for crate::alignment::record::MappingQuality {
    type Error = io::Error;

    fn try_from(raw_mapping_quality: MappingQuality<'a>) -> Result<Self, Self::Error> {
        u8::try_from(raw_mapping_quality).and_then(|n| {
            Self::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }
}
