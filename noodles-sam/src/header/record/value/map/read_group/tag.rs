//! SAM header read group tag.

use crate::header::record::value::map::tag::{self, LENGTH};

const ID: [u8; LENGTH] = [b'I', b'D'];
const BC: [u8; LENGTH] = [b'B', b'C'];
const CN: [u8; LENGTH] = [b'C', b'N'];
const DS: [u8; LENGTH] = [b'D', b'S'];
const DT: [u8; LENGTH] = [b'D', b'T'];
const FO: [u8; LENGTH] = [b'F', b'O'];
const KS: [u8; LENGTH] = [b'K', b'S'];
const LB: [u8; LENGTH] = [b'L', b'B'];
const PG: [u8; LENGTH] = [b'P', b'G'];
const PI: [u8; LENGTH] = [b'P', b'I'];
const PL: [u8; LENGTH] = [b'P', b'L'];
const PM: [u8; LENGTH] = [b'P', b'M'];
const PU: [u8; LENGTH] = [b'P', b'U'];
const SM: [u8; LENGTH] = [b'S', b'M'];

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

impl AsRef<[u8; LENGTH]> for Standard {
    fn as_ref(&self) -> &[u8; LENGTH] {
        match self {
            Standard::Id => &ID,
            Standard::Barcode => &BC,
            Standard::SequencingCenter => &CN,
            Standard::Description => &DS,
            Standard::ProducedAt => &DT,
            Standard::FlowOrder => &FO,
            Standard::KeySequence => &KS,
            Standard::Library => &LB,
            Standard::Program => &PG,
            Standard::PredictedMedianInsertSize => &PI,
            Standard::Platform => &PL,
            Standard::PlatformModel => &PM,
            Standard::PlatformUnit => &PU,
            Standard::Sample => &SM,
        }
    }
}

impl TryFrom<[u8; LENGTH]> for Standard {
    type Error = ();

    fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
        match b {
            ID => Ok(Self::Id),
            BC => Ok(Self::Barcode),
            CN => Ok(Self::SequencingCenter),
            DS => Ok(Self::Description),
            DT => Ok(Self::ProducedAt),
            FO => Ok(Self::FlowOrder),
            KS => Ok(Self::KeySequence),
            LB => Ok(Self::Library),
            PG => Ok(Self::Program),
            PI => Ok(Self::PredictedMedianInsertSize),
            PL => Ok(Self::Platform),
            PM => Ok(Self::PlatformModel),
            PU => Ok(Self::PlatformUnit),
            SM => Ok(Self::Sample),
            _ => Err(()),
        }
    }
}

impl From<Standard> for [u8; LENGTH] {
    fn from(tag: Standard) -> Self {
        match tag {
            Standard::Id => ID,
            Standard::Barcode => BC,
            Standard::SequencingCenter => CN,
            Standard::Description => DS,
            Standard::ProducedAt => DT,
            Standard::FlowOrder => FO,
            Standard::KeySequence => KS,
            Standard::Library => LB,
            Standard::Program => PG,
            Standard::PredictedMedianInsertSize => PI,
            Standard::Platform => PL,
            Standard::PlatformModel => PM,
            Standard::PlatformUnit => PU,
            Standard::Sample => SM,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_as_ref_u8_2_array_for_standard() {
        assert_eq!(Standard::Id.as_ref(), &[b'I', b'D']);
        assert_eq!(Standard::Barcode.as_ref(), &[b'B', b'C']);
        assert_eq!(Standard::SequencingCenter.as_ref(), &[b'C', b'N']);
        assert_eq!(Standard::Description.as_ref(), &[b'D', b'S']);
        assert_eq!(Standard::ProducedAt.as_ref(), &[b'D', b'T']);
        assert_eq!(Standard::FlowOrder.as_ref(), &[b'F', b'O']);
        assert_eq!(Standard::KeySequence.as_ref(), &[b'K', b'S']);
        assert_eq!(Standard::Library.as_ref(), &[b'L', b'B']);
        assert_eq!(Standard::Program.as_ref(), &[b'P', b'G']);
        assert_eq!(Standard::PredictedMedianInsertSize.as_ref(), &[b'P', b'I']);
        assert_eq!(Standard::Platform.as_ref(), &[b'P', b'L']);
        assert_eq!(Standard::PlatformModel.as_ref(), &[b'P', b'M']);
        assert_eq!(Standard::PlatformUnit.as_ref(), &[b'P', b'U']);
        assert_eq!(Standard::Sample.as_ref(), &[b'S', b'M']);
    }

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

    #[test]
    fn test_from_standard_for_u8_2_array() {
        assert_eq!(<[u8; LENGTH]>::from(Standard::Id), [b'I', b'D']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::Barcode), [b'B', b'C']);
        assert_eq!(
            <[u8; LENGTH]>::from(Standard::SequencingCenter),
            [b'C', b'N']
        );
        assert_eq!(<[u8; LENGTH]>::from(Standard::Description), [b'D', b'S']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::ProducedAt), [b'D', b'T']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::FlowOrder), [b'F', b'O']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::KeySequence), [b'K', b'S']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::Library), [b'L', b'B']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::Program), [b'P', b'G']);
        assert_eq!(
            <[u8; LENGTH]>::from(Standard::PredictedMedianInsertSize),
            [b'P', b'I']
        );
        assert_eq!(<[u8; LENGTH]>::from(Standard::Platform), [b'P', b'L']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::PlatformModel), [b'P', b'M']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::PlatformUnit), [b'P', b'U']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::Sample), [b'S', b'M']);
    }
}
