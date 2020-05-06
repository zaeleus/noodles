use std::{error, fmt, str::FromStr};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Id {
    Deletion,
    Insertion,
    Duplication,
    Inversion,
    CopyNumberVariation,
    Breakend,
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid alternative allele id: expected {{DEL, INS, DUP, INV, CNV, BND}}, got {}",
            self.0
        )
    }
}

impl FromStr for Id {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "DEL" => Ok(Self::Deletion),
            "INS" => Ok(Self::Insertion),
            "DUP" => Ok(Self::Duplication),
            "INV" => Ok(Self::Inversion),
            "CNV" => Ok(Self::CopyNumberVariation),
            "BDN" => Ok(Self::Breakend),
            _ => Err(ParseError(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("DEL".parse::<Id>()?, Id::Deletion);
        assert_eq!("INS".parse::<Id>()?, Id::Insertion);
        assert_eq!("DUP".parse::<Id>()?, Id::Duplication);
        assert_eq!("INV".parse::<Id>()?, Id::Inversion);
        assert_eq!("CNV".parse::<Id>()?, Id::CopyNumberVariation);
        assert_eq!("BDN".parse::<Id>()?, Id::Breakend);

        assert!("".parse::<Id>().is_err());
        assert!("NDL".parse::<Id>().is_err());

        Ok(())
    }
}
