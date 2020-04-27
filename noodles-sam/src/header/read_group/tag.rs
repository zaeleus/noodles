use std::{error, fmt, str::FromStr};

#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum Tag {
    Id,
    Barcode,
    SequencingCenter,
    Description,
    ProducedAt,
    FlowOrder,
    KeySequence,
    Library,
    Program,
    PredictedMedianInsertSize,
    Platform,
    PlatformModel,
    PlatformUnit,
    Sample,
    Other(String),
}

impl AsRef<str> for Tag {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => "ID",
            Self::Barcode => "BC",
            Self::SequencingCenter => "CN",
            Self::Description => "DS",
            Self::ProducedAt => "DT",
            Self::FlowOrder => "FO",
            Self::KeySequence => "KS",
            Self::Library => "LB",
            Self::Program => "PG",
            Self::PredictedMedianInsertSize => "PI",
            Self::Platform => "PL",
            Self::PlatformModel => "PM",
            Self::PlatformUnit => "PU",
            Self::Sample => "SM",
            Self::Other(s) => s,
        }
    }
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_ref())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid read group tag: '{}'", self.0)
    }
}

impl FromStr for Tag {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "ID" => Ok(Self::Id),
            "BC" => Ok(Self::Barcode),
            "CN" => Ok(Self::SequencingCenter),
            "DS" => Ok(Self::Description),
            "DT" => Ok(Self::ProducedAt),
            "FO" => Ok(Self::FlowOrder),
            "KS" => Ok(Self::KeySequence),
            "LB" => Ok(Self::Library),
            "PG" => Ok(Self::Program),
            "PI" => Ok(Self::PredictedMedianInsertSize),
            "PL" => Ok(Self::Platform),
            "PM" => Ok(Self::PlatformModel),
            "PU" => Ok(Self::PlatformUnit),
            "SM" => Ok(Self::Sample),
            _ => {
                if s.len() == 2 {
                    Ok(Self::Other(s.into()))
                } else {
                    Err(ParseError(s.into()))
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(format!("{}", Tag::Id), "ID");
        assert_eq!(format!("{}", Tag::Barcode), "BC");
        assert_eq!(format!("{}", Tag::SequencingCenter), "CN");
        assert_eq!(format!("{}", Tag::Description), "DS");
        assert_eq!(format!("{}", Tag::ProducedAt), "DT");
        assert_eq!(format!("{}", Tag::FlowOrder), "FO");
        assert_eq!(format!("{}", Tag::KeySequence), "KS");
        assert_eq!(format!("{}", Tag::Library), "LB");
        assert_eq!(format!("{}", Tag::Program), "PG");
        assert_eq!(format!("{}", Tag::PredictedMedianInsertSize), "PI");
        assert_eq!(format!("{}", Tag::Platform), "PL");
        assert_eq!(format!("{}", Tag::PlatformModel), "PM");
        assert_eq!(format!("{}", Tag::PlatformUnit), "PU");
        assert_eq!(format!("{}", Tag::Sample), "SM");
        assert_eq!(format!("{}", Tag::Other(String::from("ND"))), "ND");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("ID".parse::<Tag>()?, Tag::Id);
        assert_eq!("BC".parse::<Tag>()?, Tag::Barcode);
        assert_eq!("CN".parse::<Tag>()?, Tag::SequencingCenter);
        assert_eq!("DS".parse::<Tag>()?, Tag::Description);
        assert_eq!("DT".parse::<Tag>()?, Tag::ProducedAt);
        assert_eq!("FO".parse::<Tag>()?, Tag::FlowOrder);
        assert_eq!("KS".parse::<Tag>()?, Tag::KeySequence);
        assert_eq!("LB".parse::<Tag>()?, Tag::Library);
        assert_eq!("PG".parse::<Tag>()?, Tag::Program);
        assert_eq!("PI".parse::<Tag>()?, Tag::PredictedMedianInsertSize);
        assert_eq!("PL".parse::<Tag>()?, Tag::Platform);
        assert_eq!("PM".parse::<Tag>()?, Tag::PlatformModel);
        assert_eq!("PU".parse::<Tag>()?, Tag::PlatformUnit);
        assert_eq!("SM".parse::<Tag>()?, Tag::Sample);

        assert_eq!("ND".parse::<Tag>()?, Tag::Other(String::from("ND")));

        assert!("".parse::<Tag>().is_err());
        assert!("NDL".parse::<Tag>().is_err());

        Ok(())
    }
}
