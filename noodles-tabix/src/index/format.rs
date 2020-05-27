use std::{convert::TryFrom, error, fmt};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    Generic,
    Sam,
    Vcf,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TryFromIntError(i32);

impl error::Error for TryFromIntError {}

impl fmt::Display for TryFromIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid format: expected 0..=2, got {}", self.0)
    }
}

impl TryFrom<i32> for Format {
    type Error = TryFromIntError;

    fn try_from(n: i32) -> Result<Self, Self::Error> {
        match n {
            0 => Ok(Self::Generic),
            1 => Ok(Self::Sam),
            2 => Ok(Self::Vcf),
            _ => Err(TryFromIntError(n)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_i32_for_format() {
        assert_eq!(Format::try_from(0), Ok(Format::Generic));
        assert_eq!(Format::try_from(1), Ok(Format::Sam));
        assert_eq!(Format::try_from(2), Ok(Format::Vcf));
        assert_eq!(Format::try_from(3), Err(TryFromIntError(3)));
    }
}
