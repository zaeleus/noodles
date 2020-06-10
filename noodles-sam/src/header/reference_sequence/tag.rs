use std::{error, fmt, str::FromStr};

/// A SAM header reference sequence tag.
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum Tag {
    Name,
    Length,
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

impl AsRef<str> for Tag {
    fn as_ref(&self) -> &str {
        match self {
            Self::Name => "SN",
            Self::Length => "LN",
            Self::AlternativeLocus => "AH",
            Self::AlternativeNames => "AN",
            Self::AssemblyId => "AS",
            Self::Description => "DS",
            Self::Md5Checksum => "M5",
            Self::Species => "SP",
            Self::MoleculeTopology => "TP",
            Self::Uri => "UR",
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
        write!(f, "invalid reference sequence tag: '{}'", self.0)
    }
}

impl FromStr for Tag {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "SN" => Ok(Self::Name),
            "LN" => Ok(Self::Length),
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
        assert_eq!(format!("{}", Tag::Name), "SN");
        assert_eq!(format!("{}", Tag::Length), "LN");
        assert_eq!(format!("{}", Tag::AlternativeLocus), "AH");
        assert_eq!(format!("{}", Tag::AlternativeNames), "AN");
        assert_eq!(format!("{}", Tag::AssemblyId), "AS");
        assert_eq!(format!("{}", Tag::Description), "DS");
        assert_eq!(format!("{}", Tag::Md5Checksum), "M5");
        assert_eq!(format!("{}", Tag::Species), "SP");
        assert_eq!(format!("{}", Tag::MoleculeTopology), "TP");
        assert_eq!(format!("{}", Tag::Uri), "UR");
        assert_eq!(format!("{}", Tag::Other(String::from("ND"))), "ND");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("SN".parse::<Tag>()?, Tag::Name);
        assert_eq!("LN".parse::<Tag>()?, Tag::Length);
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
