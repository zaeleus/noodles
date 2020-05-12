mod base;

pub use self::base::Base;

use std::{convert::TryFrom, error, fmt, ops::Deref, str::FromStr};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReferenceBases(Vec<Base>);

impl Deref for ReferenceBases {
    type Target = [Base];

    fn deref(&self) -> &Self::Target {
        &self.0
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
        s.chars()
            .map(|c| c.to_ascii_uppercase())
            .map(Base::try_from)
            .collect::<Result<_, _>>()
            .map(ReferenceBases)
            .map_err(|_| ParseError(s.into()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let expected = [Base::A, Base::T, Base::C, Base::G, Base::N];

        let bases: ReferenceBases = "ATCGN".parse()?;
        assert_eq!(&bases[..], &expected[..]);

        let bases: ReferenceBases = "atcgn".parse()?;
        assert_eq!(&bases[..], &expected[..]);

        let bases: ReferenceBases = "AtCgN".parse()?;
        assert_eq!(&bases[..], &expected[..]);

        Ok(())
    }
}
