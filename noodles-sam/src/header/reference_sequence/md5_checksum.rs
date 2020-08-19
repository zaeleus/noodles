use std::{error, fmt, num, ops::Deref, str::FromStr};

/// A SAM header reference sequence MD5 checksum.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Md5Checksum([u8; 16]);

impl Deref for Md5Checksum {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Md5Checksum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for byte in self.iter() {
            write!(f, "{:x}", byte)?;
        }

        Ok(())
    }
}

impl From<[u8; 16]> for Md5Checksum {
    fn from(checksum: [u8; 16]) -> Self {
        Self(checksum)
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidLength(usize),
    InvalidHex(num::ParseIntError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLength(len) => write!(f, "expected length to be 32, got {}", len),
            Self::InvalidHex(e) => write!(f, "invalid hex: {}", e),
        }
    }
}

impl FromStr for Md5Checksum {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() != 32 {
            return Err(ParseError::InvalidLength(s.len()));
        }

        let mut checksum = [0; 16];

        for (i, value) in checksum.iter_mut().enumerate() {
            let start = 2 * i;
            let end = start + 1;
            let hex = &s[start..=end];
            *value = u8::from_str_radix(hex, 16).map_err(ParseError::InvalidHex)?;
        }

        Ok(Self(checksum))
    }
}

impl From<Md5Checksum> for [u8; 16] {
    fn from(md5_checksum: Md5Checksum) -> Self {
        md5_checksum.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let md5_checksum = Md5Checksum::from([
            0xd7, 0xeb, 0xa3, 0x11, 0x42, 0x1b, 0xbc, 0x9d, 0x3a, 0xda, 0x44, 0x70, 0x9d, 0xd6,
            0x15, 0x34,
        ]);

        assert_eq!(md5_checksum.to_string(), "d7eba311421bbc9d3ada44709dd61534");
    }

    #[test]
    fn test_from_u8_array_for_md5_checksum() {
        let checksum = [
            0xd7, 0xeb, 0xa3, 0x11, 0x42, 0x1b, 0xbc, 0x9d, 0x3a, 0xda, 0x44, 0x70, 0x9d, 0xd6,
            0x15, 0x34,
        ];
        assert_eq!(*Md5Checksum::from(checksum), checksum);
    }

    #[test]
    fn test_parse() {
        assert_eq!(
            "d7eba311421bbc9d3ada44709dd61534".parse(),
            Ok(Md5Checksum::from([
                0xd7, 0xeb, 0xa3, 0x11, 0x42, 0x1b, 0xbc, 0x9d, 0x3a, 0xda, 0x44, 0x70, 0x9d, 0xd6,
                0x15, 0x34,
            ]))
        );

        assert_eq!("".parse::<Md5Checksum>(), Err(ParseError::InvalidLength(0)));
        assert_eq!(
            "d7".parse::<Md5Checksum>(),
            Err(ParseError::InvalidLength(2))
        );
        assert_eq!(
            "838f8d9acc45bd36e3213c47c3222e644f44c959fa370bbfa6df46b171c02f0c"
                .parse::<Md5Checksum>(),
            Err(ParseError::InvalidLength(64))
        );

        assert!(matches!(
            "n7eba311421bbc9d3ada44709dd61534".parse::<Md5Checksum>(),
            Err(ParseError::InvalidHex(_))
        ));
    }

    #[test]
    fn test_from_md5_checksum_for_u8_array() {
        let checksum = [
            0xd7, 0xeb, 0xa3, 0x11, 0x42, 0x1b, 0xbc, 0x9d, 0x3a, 0xda, 0x44, 0x70, 0x9d, 0xd6,
            0x15, 0x34,
        ];
        let md5_checksum = Md5Checksum::from(checksum);
        assert_eq!(<[u8; 16]>::from(md5_checksum), checksum);
    }
}
