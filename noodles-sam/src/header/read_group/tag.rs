use std::str::FromStr;

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

impl FromStr for Tag {
    type Err = ();

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
                    Err(())
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ()> {
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
