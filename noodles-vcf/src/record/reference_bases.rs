mod base;

pub use self::base::Base;

use std::{convert::TryFrom, error, fmt, ops::Deref, str::FromStr};

use super::MISSING_FIELD;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReferenceBases(Vec<Base>);

impl Deref for ReferenceBases {
    type Target = [Base];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for ReferenceBases {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for base in self.iter() {
            write!(f, "{}", char::from(*base))?;
        }

        Ok(())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid reference bases: {}", self.0)
    }
}

impl FromStr for ReferenceBases {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" | MISSING_FIELD => Err(ParseError(s.into())),
            _ => s
                .chars()
                .map(|c| c.to_ascii_uppercase())
                .map(Base::try_from)
                .collect::<Result<_, _>>()
                .map(ReferenceBases)
                .map_err(|_| ParseError(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let reference_bases = ReferenceBases(vec![Base::A, Base::T, Base::C, Base::G, Base::N]);
        assert_eq!(reference_bases.to_string(), "ATCGN");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let expected = [Base::A, Base::T, Base::C, Base::G, Base::N];

        let bases: ReferenceBases = "ATCGN".parse()?;
        assert_eq!(&bases[..], &expected[..]);

        let bases: ReferenceBases = "atcgn".parse()?;
        assert_eq!(&bases[..], &expected[..]);

        let bases: ReferenceBases = "AtCgN".parse()?;
        assert_eq!(&bases[..], &expected[..]);

        assert!("".parse::<ReferenceBases>().is_err());
        assert!(".".parse::<ReferenceBases>().is_err());

        Ok(())
    }
}
