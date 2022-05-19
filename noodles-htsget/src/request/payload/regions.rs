use noodles_core::Region;
use serde::{
    ser::{SerializeMap, SerializeSeq},
    Serialize, Serializer,
};

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Regions(pub(super) Vec<Region>);

impl Regions {
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl Serialize for Regions {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(self.0.len()))?;

        for region in &self.0 {
            seq.serialize_element(&HtsgetRegion(region))?;
        }

        seq.end()
    }
}

struct HtsgetRegion<'a>(&'a Region);

impl<'a> Serialize for HtsgetRegion<'a> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use crate::resolve_interval;

        let mut map = serializer.serialize_map(None)?;

        map.serialize_entry("referenceName", self.0.name())?;

        let (resolved_start, resolved_end) = resolve_interval(self.0.interval());

        if let Some(position) = resolved_start {
            let start = u32::try_from(position).unwrap();
            map.serialize_entry("start", &start)?;
        }

        if let Some(position) = resolved_end {
            let end = u32::try_from(position).unwrap();
            map.serialize_entry("end", &end)?;
        }

        map.end()
    }
}
