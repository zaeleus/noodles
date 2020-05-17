use std::{error, fmt, str::FromStr};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Type {
    Deletion,
    Insertion,
    Duplication,
    Inversion,
    CopyNumberVariation,
    Breakend,
}

impl AsRef<str> for Type {
    fn as_ref(&self) -> &str {
        match self {
            Self::Deletion => "DEL",
            Self::Insertion => "INS",
            Self::Duplication => "DUP",
            Self::Inversion => "INV",
            Self::CopyNumberVariation => "CNV",
            Self::Breakend => "BND",
        }
    }
}

impl fmt::Display for Type {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid alternate bases symbol type: expected {{DEL, INS, DUP, INV, CNV, BND}}, got {}",
            self.0
        )
    }
}

impl FromStr for Type {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "DEL" => Ok(Self::Deletion),
            "INS" => Ok(Self::Insertion),
            "DUP" => Ok(Self::Duplication),
            "INV" => Ok(Self::Inversion),
            "CNV" => Ok(Self::CopyNumberVariation),
            "BND" => Ok(Self::Breakend),
            _ => Err(ParseError(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Type::Deletion.to_string(), "DEL");
        assert_eq!(Type::Insertion.to_string(), "INS");
        assert_eq!(Type::Duplication.to_string(), "DUP");
        assert_eq!(Type::Inversion.to_string(), "INV");
        assert_eq!(Type::CopyNumberVariation.to_string(), "CNV");
        assert_eq!(Type::Breakend.to_string(), "BND");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("DEL".parse::<Type>()?, Type::Deletion);
        assert_eq!("INS".parse::<Type>()?, Type::Insertion);
        assert_eq!("DUP".parse::<Type>()?, Type::Duplication);
        assert_eq!("INV".parse::<Type>()?, Type::Inversion);
        assert_eq!("CNV".parse::<Type>()?, Type::CopyNumberVariation);
        assert_eq!("BND".parse::<Type>()?, Type::Breakend);

        assert!("".parse::<Type>().is_err());
        assert!("NDL".parse::<Type>().is_err());

        Ok(())
    }
}
