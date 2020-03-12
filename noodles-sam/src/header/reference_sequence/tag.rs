use std::str::FromStr;

#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum Tag {
    Name,
    Len,
    AlternativeLocus,
    AlternativeNames,
    AssemblyId,
    Description,
    Md5Checksum,
    Species,
    MoleculeTopology,
    Uri,
    Other(String),
}

impl FromStr for Tag {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "SN" => Ok(Self::Name),
            "LN" => Ok(Self::Len),
            "AH" => Ok(Self::AlternativeLocus),
            "AN" => Ok(Self::AlternativeNames),
            "AS" => Ok(Self::AssemblyId),
            "DS" => Ok(Self::Description),
            "M5" => Ok(Self::Md5Checksum),
            "SP" => Ok(Self::Species),
            "TP" => Ok(Self::MoleculeTopology),
            "UR" => Ok(Self::Uri),
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
        assert_eq!("SN".parse::<Tag>()?, Tag::Name);
        assert_eq!("LN".parse::<Tag>()?, Tag::Len);
        assert_eq!("AH".parse::<Tag>()?, Tag::AlternativeLocus);
        assert_eq!("AN".parse::<Tag>()?, Tag::AlternativeNames);
        assert_eq!("AS".parse::<Tag>()?, Tag::AssemblyId);
        assert_eq!("DS".parse::<Tag>()?, Tag::Description);
        assert_eq!("M5".parse::<Tag>()?, Tag::Md5Checksum);
        assert_eq!("SP".parse::<Tag>()?, Tag::Species);
        assert_eq!("TP".parse::<Tag>()?, Tag::MoleculeTopology);
        assert_eq!("UR".parse::<Tag>()?, Tag::Uri);

        assert_eq!("ND".parse::<Tag>()?, Tag::Other(String::from("ND")));

        assert!("".parse::<Tag>().is_err());
        assert!("NDL".parse::<Tag>().is_err());

        Ok(())
    }
}
