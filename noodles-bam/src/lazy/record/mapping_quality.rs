/// A raw BAM record mapping quality.
#[derive(Debug, Eq, PartialEq)]
pub struct MappingQuality(u8);

impl MappingQuality {
    pub(super) fn new(n: u8) -> Self {
        Self(n)
    }
}

impl From<MappingQuality> for u8 {
    fn from(mapping_quality: MappingQuality) -> Self {
        mapping_quality.0
    }
}
