use std::{fmt, str::FromStr};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
    Irrelevant,
    Unknown,
}

impl Default for Strand {
    fn default() -> Self {
        Self::Irrelevant
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = match self {
            Self::Forward => '+',
            Self::Reverse => '-',
            Self::Irrelevant => '.',
            Self::Unknown => '?',
        };

        write!(f, "{}", s)
    }
}

impl FromStr for Strand {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Self::Forward),
            "-" => Ok(Self::Reverse),
            "." => Ok(Self::Irrelevant),
            "?" => Ok(Self::Unknown),
            _ => Err(format!("invalid strand '{}'", s)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Strand;

    #[test]
    fn test_default() {
        assert_eq!(Strand::default(), Strand::Irrelevant);
    }

    #[test]
    fn test_display() {
        assert_eq!(format!("{}", Strand::Forward), "+");
        assert_eq!(format!("{}", Strand::Reverse), "-");
        assert_eq!(format!("{}", Strand::Irrelevant), ".");
        assert_eq!(format!("{}", Strand::Unknown), "?");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("+".parse(), Ok(Strand::Forward));
        assert_eq!("-".parse(), Ok(Strand::Reverse));
        assert_eq!(".".parse(), Ok(Strand::Irrelevant));
        assert_eq!("?".parse(), Ok(Strand::Unknown));
        assert!("!".parse::<Strand>().is_err());
    }
}
