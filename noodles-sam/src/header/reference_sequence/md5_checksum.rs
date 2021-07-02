//! SAM header reference sequence MD5 checksum.

use std::{error, fmt, ops::Deref, str::FromStr};

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
            write!(f, "{:02x}", byte)?;
        }

        Ok(())
    }
}

impl From<[u8; 16]> for Md5Checksum {
    fn from(checksum: [u8; 16]) -> Self {
        Self(checksum)
    }
}

/// An error returned when a raw SAM header reference sequence MD5 checksum fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The length is invalid.
    InvalidLength(usize),
    /// The input has an invalid hex digit.
    InvalidHexDigit(char),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLength(len) => write!(f, "expected length to be 32, got {}", len),
            Self::InvalidHexDigit(c) => write!(f, "invalid hex digit: {}", c),
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

        for (digits, value) in s.as_bytes().chunks(2).zip(checksum.iter_mut()) {
            let l = parse_digit(digits[0])?;
            let r = parse_digit(digits[1])?;
            *value = l << 4 | r;
        }

        Ok(Self(checksum))
    }
}

impl From<Md5Checksum> for [u8; 16] {
    fn from(md5_checksum: Md5Checksum) -> Self {
        md5_checksum.0
    }
}

// § 1.3.2 Reference MD5 calculation (2021-06-03): "The MD5 digest is ... presented as a 32
// character lowercase hexadecimal number."
fn parse_digit(b: u8) -> Result<u8, ParseError> {
    match b {
        b'a'..=b'f' => Ok(b - b'a' + 10),
        b'0'..=b'9' => Ok(b - b'0'),
        _ => Err(ParseError::InvalidHexDigit(char::from(b))),
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

        let md5_checksum = Md5Checksum::from([
            0xb0, 0x0c, 0x61, 0xdf, 0xed, 0x4a, 0x92, 0xfd, 0xfb, 0x24, 0x4d, 0x35, 0x79, 0x05,
            0x56, 0xeb,
        ]);
        assert_eq!(md5_checksum.to_string(), "b00c61dfed4a92fdfb244d35790556eb");
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

        assert_eq!(
            "D7EBA311421BBC9D3ADA44709DD61534".parse::<Md5Checksum>(),
            Err(ParseError::InvalidHexDigit('D'))
        );

        assert!(matches!(
            "n7eba311421bbc9d3ada44709dd61534".parse::<Md5Checksum>(),
            Err(ParseError::InvalidHexDigit('n'))
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
