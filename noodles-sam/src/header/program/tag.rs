use std::str::FromStr;

#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum Tag {
    Id,
    Name,
    CommandLine,
    PreviousId,
    Description,
    Version,
    Other(String),
}

impl FromStr for Tag {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "ID" => Ok(Self::Id),
            "PN" => Ok(Self::Name),
            "CL" => Ok(Self::CommandLine),
            "PP" => Ok(Self::PreviousId),
            "DS" => Ok(Self::Description),
            "VN" => Ok(Self::Version),
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
        assert_eq!("PN".parse::<Tag>()?, Tag::Name);
        assert_eq!("CL".parse::<Tag>()?, Tag::CommandLine);
        assert_eq!("PP".parse::<Tag>()?, Tag::PreviousId);
        assert_eq!("DS".parse::<Tag>()?, Tag::Description);
        assert_eq!("VN".parse::<Tag>()?, Tag::Version);

        assert_eq!("ND".parse::<Tag>()?, Tag::Other(String::from("ND")));

        assert!("".parse::<Tag>().is_err());
        assert!("NDL".parse::<Tag>().is_err());

        Ok(())
    }
}
