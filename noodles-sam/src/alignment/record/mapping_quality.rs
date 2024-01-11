use std::io;

/// An alignment record mapping quality.
pub trait MappingQuality {
    /// Converts a mapping quality to a `u8`.
    fn try_to_u8(&self) -> io::Result<u8>;
}

impl TryFrom<&dyn MappingQuality> for u8 {
    type Error = io::Error;

    fn try_from(raw_mapping_quality: &dyn MappingQuality) -> Result<Self, Self::Error> {
        raw_mapping_quality.try_to_u8()
    }
}

impl TryFrom<&dyn MappingQuality> for crate::alignment::record_buf::MappingQuality {
    type Error = io::Error;

    fn try_from(raw_mapping_quality: &dyn MappingQuality) -> Result<Self, Self::Error> {
        u8::try_from(raw_mapping_quality).and_then(|n| {
            Self::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }
}
