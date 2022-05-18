use noodles_core::Region;

use super::Kind;
use crate::Format;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Payload {
    format: Format,
    regions: Vec<Region>,
}

impl Payload {
    pub fn regions(&self) -> &[Region] {
        &self.regions
    }

    pub fn regions_mut(&mut self) -> &mut Vec<Region> {
        &mut self.regions
    }
}

impl From<Kind> for Payload {
    fn from(kind: Kind) -> Self {
        let format = match kind {
            Kind::Reads => Format::Bam,
            Kind::Variants => Format::Vcf,
        };

        Self {
            format,
            regions: Vec::new(),
        }
    }
}
