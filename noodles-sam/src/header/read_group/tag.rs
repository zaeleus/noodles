//! SAM header read group tag.

use std::{
    convert::TryFrom,
    error,
    fmt::{self, Write},
    str::FromStr,
};

const LENGTH: usize = 2;

/// A SAM header read group tag.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum Tag {
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
    /// Any other read group tag.
    Other([u8; LENGTH]),
}

impl AsRef<[u8; LENGTH]> for Tag {
    fn as_ref(&self) -> &[u8; LENGTH] {
        match self {
            Self::Id => b"ID",
            Self::Barcode => b"BC",
            Self::SequencingCenter => b"CN",
            Self::Description => b"DS",
            Self::ProducedAt => b"DT",
            Self::FlowOrder => b"FO",
            Self::KeySequence => b"KS",
            Self::Library => b"LB",
            Self::Program => b"PG",
            Self::PredictedMedianInsertSize => b"PI",
            Self::Platform => b"PL",
            Self::PlatformModel => b"PM",
            Self::PlatformUnit => b"PU",
            Self::Sample => b"SM",
            Self::Other(tag) => tag,
        }
    }
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let bytes = self.as_ref();
        f.write_char(char::from(bytes[0]))?;
        f.write_char(char::from(bytes[1]))?;
        Ok(())
    }
}

/// An error returned when a raw SAM header read group tag fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Tag {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let bytes = s.as_bytes();

        match bytes.len() {
            0 => Err(ParseError::Empty),
            2 => Self::try_from([bytes[0], bytes[1]]),
            _ => Err(ParseError::Invalid),
        }
    }
}

impl TryFrom<[u8; LENGTH]> for Tag {
    type Error = ParseError;

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
            _ => Ok(Self::Other(b)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Tag::Id.to_string(), "ID");
        assert_eq!(Tag::Barcode.to_string(), "BC");
        assert_eq!(Tag::SequencingCenter.to_string(), "CN");
        assert_eq!(Tag::Description.to_string(), "DS");
        assert_eq!(Tag::ProducedAt.to_string(), "DT");
        assert_eq!(Tag::FlowOrder.to_string(), "FO");
        assert_eq!(Tag::KeySequence.to_string(), "KS");
        assert_eq!(Tag::Library.to_string(), "LB");
        assert_eq!(Tag::Program.to_string(), "PG");
        assert_eq!(Tag::PredictedMedianInsertSize.to_string(), "PI");
        assert_eq!(Tag::Platform.to_string(), "PL");
        assert_eq!(Tag::PlatformModel.to_string(), "PM");
        assert_eq!(Tag::PlatformUnit.to_string(), "PU");
        assert_eq!(Tag::Sample.to_string(), "SM");
        assert_eq!(Tag::Other([b'N', b'D']).to_string(), "ND");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("ID".parse(), Ok(Tag::Id));
        assert_eq!("BC".parse(), Ok(Tag::Barcode));
        assert_eq!("CN".parse(), Ok(Tag::SequencingCenter));
        assert_eq!("DS".parse(), Ok(Tag::Description));
        assert_eq!("DT".parse(), Ok(Tag::ProducedAt));
        assert_eq!("FO".parse(), Ok(Tag::FlowOrder));
        assert_eq!("KS".parse(), Ok(Tag::KeySequence));
        assert_eq!("LB".parse(), Ok(Tag::Library));
        assert_eq!("PG".parse(), Ok(Tag::Program));
        assert_eq!("PI".parse(), Ok(Tag::PredictedMedianInsertSize));
        assert_eq!("PL".parse(), Ok(Tag::Platform));
        assert_eq!("PM".parse(), Ok(Tag::PlatformModel));
        assert_eq!("PU".parse(), Ok(Tag::PlatformUnit));
        assert_eq!("SM".parse(), Ok(Tag::Sample));
        assert_eq!("ND".parse(), Ok(Tag::Other([b'N', b'D'])));

        assert_eq!("".parse::<Tag>(), Err(ParseError::Empty));
        assert_eq!("NDL".parse::<Tag>(), Err(ParseError::Invalid));
    }
}
