//! SAM header read group tag.

use crate::header::record::value::map::tag::{self, LENGTH};

/// A SAM header read group tag.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum Standard {
    /// Read group ID (`ID`).
    Id,
    /// Barcode sequence (`BC`).
    Barcode,
    /// Sequencing center (`CN`).
    SequencingCenter,
    /// Description (`DS`).
    Description,
    /// Datetime of run (`DT`).
    ProducedAt,
    /// Flow order (`FO`).
    FlowOrder,
    /// Key sequence (`KS`).
    KeySequence,
    /// Library (`LB`).
    Library,
    /// Programs used (`PG`).
    Program,
    /// Predicted median insert size (`PI`).
    PredictedMedianInsertSize,
    /// Platform used (`PL`).
    Platform,
    /// Platform model (`PM`).
    PlatformModel,
    /// Platform unit (`PU`).
    PlatformUnit,
    /// Sample (`SM`).
    Sample,
}

impl tag::Standard for Standard {}

impl TryFrom<[u8; LENGTH]> for Standard {
    type Error = ();

    fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
        match &b {
            b"ID" => Ok(Self::Id),
            b"BC" => Ok(Self::Barcode),
            b"CN" => Ok(Self::SequencingCenter),
            b"DS" => Ok(Self::Description),
            b"DT" => Ok(Self::ProducedAt),
            b"FO" => Ok(Self::FlowOrder),
            b"KS" => Ok(Self::KeySequence),
            b"LB" => Ok(Self::Library),
            b"PG" => Ok(Self::Program),
            b"PI" => Ok(Self::PredictedMedianInsertSize),
            b"PL" => Ok(Self::Platform),
            b"PM" => Ok(Self::PlatformModel),
            b"PU" => Ok(Self::PlatformUnit),
            b"SM" => Ok(Self::Sample),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_u8_array_for_standard() {
        assert_eq!(Standard::try_from([b'I', b'D']), Ok(Standard::Id));
        assert_eq!(Standard::try_from([b'B', b'C']), Ok(Standard::Barcode));
        assert_eq!(
            Standard::try_from([b'C', b'N']),
            Ok(Standard::SequencingCenter)
        );
        assert_eq!(Standard::try_from([b'D', b'S']), Ok(Standard::Description));
        assert_eq!(Standard::try_from([b'D', b'T']), Ok(Standard::ProducedAt));
        assert_eq!(Standard::try_from([b'F', b'O']), Ok(Standard::FlowOrder));
        assert_eq!(Standard::try_from([b'K', b'S']), Ok(Standard::KeySequence));
        assert_eq!(Standard::try_from([b'L', b'B']), Ok(Standard::Library));
        assert_eq!(Standard::try_from([b'P', b'G']), Ok(Standard::Program));
        assert_eq!(
            Standard::try_from([b'P', b'I']),
            Ok(Standard::PredictedMedianInsertSize)
        );
        assert_eq!(Standard::try_from([b'P', b'L']), Ok(Standard::Platform));
        assert_eq!(
            Standard::try_from([b'P', b'M']),
            Ok(Standard::PlatformModel)
        );
        assert_eq!(Standard::try_from([b'P', b'U']), Ok(Standard::PlatformUnit));
        assert_eq!(Standard::try_from([b'S', b'M']), Ok(Standard::Sample));

        assert_eq!(Standard::try_from([b'N', b'D']), Err(()));
    }
}
