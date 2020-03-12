use std::str::FromStr;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    Header,
    ReferenceSequence,
    ReadGroup,
    Program,
    Comment,
}

impl FromStr for Kind {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "@HD" => Ok(Self::Header),
            "@SQ" => Ok(Self::ReferenceSequence),
            "@RG" => Ok(Self::ReadGroup),
            "@PG" => Ok(Self::Program),
            "@CO" => Ok(Self::Comment),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ()> {
        assert_eq!("@HD".parse::<Kind>()?, Kind::Header);
        assert_eq!("@SQ".parse::<Kind>()?, Kind::ReferenceSequence);
        assert_eq!("@RG".parse::<Kind>()?, Kind::ReadGroup);
        assert_eq!("@PG".parse::<Kind>()?, Kind::Program);
        assert_eq!("@CO".parse::<Kind>()?, Kind::Comment);

        assert!("@NO".parse::<Kind>().is_err());
        assert!("HD".parse::<Kind>().is_err());

        Ok(())
    }
}
