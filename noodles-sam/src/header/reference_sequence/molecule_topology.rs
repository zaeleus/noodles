//! SAM header reference sequence molecule topology.

use std::{error, fmt, str::FromStr};

/// A SAM header reference sequence molecule topology (`TP`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MoleculeTopology {
    /// Linear molecule topology (`linear`).
    Linear,
    /// Circular/cyclic molecule topology (`circular`).
    Circular,
}

impl AsRef<str> for MoleculeTopology {
    fn as_ref(&self) -> &str {
        match self {
            Self::Linear => "linear",
            Self::Circular => "circular",
        }
    }
}

impl Default for MoleculeTopology {
    fn default() -> Self {
        Self::Linear
    }
}

impl fmt::Display for MoleculeTopology {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw SAM header reference sequence molecule topology fails to parse.
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

impl FromStr for MoleculeTopology {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            "linear" => Ok(Self::Linear),
            "circular" => Ok(Self::Circular),
            _ => Err(ParseError::Invalid),
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
    fn test_fmt() {
        assert_eq!(MoleculeTopology::Linear.to_string(), "linear");
        assert_eq!(MoleculeTopology::Circular.to_string(), "circular");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("linear".parse(), Ok(MoleculeTopology::Linear));
        assert_eq!("circular".parse(), Ok(MoleculeTopology::Circular));

        assert_eq!("".parse::<MoleculeTopology>(), Err(ParseError::Empty));
        assert_eq!(
            "noodles".parse::<MoleculeTopology>(),
            Err(ParseError::Invalid)
        );
        assert_eq!(
            "Linear".parse::<MoleculeTopology>(),
            Err(ParseError::Invalid)
        );
    }
}
