use std::{error, fmt, str::FromStr};

/// A SAM header reference sequence molecule topology (`TP`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MoleculeTopology {
    /// Linear molecule topology (`linear`).
    Linear,
    /// Circular/cyclic molecule topology (`circular`).
    Circular,
}

impl Default for MoleculeTopology {
    fn default() -> Self {
        Self::Linear
    }
}

/// An error returned when a raw SAM header reference sequence molecule topology fails to parse.
#[derive(Debug, Eq, PartialEq)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid molecule topology: expected {{linear, circular}}, got {}",
            self.0
        )
    }
}

impl FromStr for MoleculeTopology {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "linear" => Ok(Self::Linear),
            "circular" => Ok(Self::Circular),
            _ => Err(ParseError(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(MoleculeTopology::default(), MoleculeTopology::Linear);
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("linear".parse(), Ok(MoleculeTopology::Linear));
        assert_eq!("circular".parse(), Ok(MoleculeTopology::Circular));

        assert_eq!(
            "".parse::<MoleculeTopology>(),
            Err(ParseError(String::from("")))
        );

        assert_eq!(
            "noodles".parse::<MoleculeTopology>(),
            Err(ParseError(String::from("noodles")))
        );

        assert_eq!(
            "Linear".parse::<MoleculeTopology>(),
            Err(ParseError(String::from("Linear")))
        );

        Ok(())
    }
}
