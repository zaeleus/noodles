use std::str;

use noodles_core::{region::Interval, Region};
use serde::{
    ser::{self, SerializeMap, SerializeSeq},
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
        let mut map = serializer.serialize_map(None)?;

        let name = str::from_utf8(self.0.name()).map_err(ser::Error::custom)?;
        map.serialize_entry("referenceName", name)?;

        let (resolved_start, resolved_end) = resolve_interval(self.0.interval());

        if let Some(position) = resolved_start {
            let start = u32::try_from(position).map_err(ser::Error::custom)?;
            map.serialize_entry("start", &start)?;
        }

        if let Some(position) = resolved_end {
            let end = u32::try_from(position).map_err(ser::Error::custom)?;
            map.serialize_entry("end", &end)?;
        }

        map.end()
    }
}

fn resolve_interval<I>(interval: I) -> (Option<usize>, Option<usize>)
where
    I: Into<Interval>,
{
    let interval = interval.into();
    let start = interval.start().map(|position| usize::from(position) - 1);
    let end = interval.end().map(usize::from);
    (start, end)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resolve_interval() -> std::result::Result<(), noodles_core::position::TryFromIntError> {
        use noodles_core::Position;

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;

        assert_eq!(resolve_interval(start..=end), (Some(7), Some(13)));
        assert_eq!(resolve_interval(start..), (Some(7), None));
        assert_eq!(resolve_interval(..=end), (None, Some(13)));
        assert_eq!(resolve_interval(..), (None, None));

        Ok(())
    }
}
