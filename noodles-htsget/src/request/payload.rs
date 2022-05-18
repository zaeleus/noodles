mod regions;

use noodles_core::Region;
use serde::Serialize;

use self::regions::Regions;
use super::Kind;
use crate::Format;

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Payload {
    format: Format,
    regions: Regions,
}

impl Payload {
    pub fn regions_mut(&mut self) -> &mut Vec<Region> {
        &mut self.regions.0
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
            regions: Regions::default(),
        }
    }
}
