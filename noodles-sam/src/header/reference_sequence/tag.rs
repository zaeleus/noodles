use std::{error, fmt, str::FromStr};

/// A SAM header reference sequence tag.
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum Tag {
    /// Reference sequence name (`SN`).
    Name,
    /// Reference sequence length (`LN`).
    Length,
    /// Alternate locus (`AH`).
    AlternativeLocus,
    /// Alternate reference sequence names (`AN`).
    AlternativeNames,
    /// Genome assembly ID (`AS`).
    AssemblyId,
    /// Description (`DS`).
    Description,
    /// MD5 checksum of the reference sequence (`M5`).
    Md5Checksum,
    /// Species (`SP`).
    Species,
    /// Molecule topology (`TP`).
    MoleculeTopology,
    /// URI of the reference sequence (`UR`).
    Uri,
    /// Any other reference sequence tag.
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

/// An error returned when a raw SAM header reference sequence tag fails to parse.
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
        match s {
            "" => Err(ParseError::Empty),
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
                    Err(ParseError::Invalid)
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
    fn test_from_str() {
        assert_eq!("SN".parse::<Tag>(), Ok(Tag::Name));
        assert_eq!("LN".parse::<Tag>(), Ok(Tag::Length));
        assert_eq!("AH".parse::<Tag>(), Ok(Tag::AlternativeLocus));
        assert_eq!("AN".parse::<Tag>(), Ok(Tag::AlternativeNames));
        assert_eq!("AS".parse::<Tag>(), Ok(Tag::AssemblyId));
        assert_eq!("DS".parse::<Tag>(), Ok(Tag::Description));
        assert_eq!("M5".parse::<Tag>(), Ok(Tag::Md5Checksum));
        assert_eq!("SP".parse::<Tag>(), Ok(Tag::Species));
        assert_eq!("TP".parse::<Tag>(), Ok(Tag::MoleculeTopology));
        assert_eq!("UR".parse::<Tag>(), Ok(Tag::Uri));
        assert_eq!("ND".parse::<Tag>(), Ok(Tag::Other(String::from("ND"))));

        assert_eq!("".parse::<Tag>(), Err(ParseError::Empty));
        assert_eq!("NDL".parse::<Tag>(), Err(ParseError::Invalid));
    }
}
