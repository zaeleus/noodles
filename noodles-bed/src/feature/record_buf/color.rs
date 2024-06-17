//! BED record color.

use std::{error, fmt, num, str::FromStr};

const DELIMITER: char = ',';

/// A BED record color.
///
/// A color is represented as an RGB triplet, where each component ranges from 0 to 255, inclusive.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Color {
    r: u8,
    g: u8,
    b: u8,
}

impl Color {
    /// White `(255, 255, 255)`.
    pub const WHITE: Self = Self::new(0xff, 0xff, 0xff);

    /// Black `(0, 0, 0)`.
    pub const BLACK: Self = Self::new(0x00, 0x00, 0x00);

    /// Red `(255, 0, 0)`.
    pub const RED: Self = Self::new(0xff, 0x00, 0x00);

    /// Lime `(0, 255, 0)`.
    pub const LIME: Self = Self::new(0x00, 0xff, 0x00);

    /// Green `(0, 128, 0)`.
    pub const GREEN: Self = Self::new(0x00, 0x80, 0x00);

    /// Blue `(0, 0, 255)`.
    pub const BLUE: Self = Self::new(0x00, 0x00, 0xff);

    /// Creates a BED record color.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::feature::record_buf::Color;
    /// let color = Color::new(250, 128, 114);
    /// ```
    pub const fn new(r: u8, g: u8, b: u8) -> Self {
        Self { r, g, b }
    }

    /// Returns the value of the red component.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::feature::record_buf::Color;
    /// let color = Color::new(250, 128, 114);
    /// assert_eq!(color.red(), 250);
    /// ```
    pub const fn red(&self) -> u8 {
        self.r
    }

    /// Returns the value of the red component.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::feature::record_buf::Color;
    /// let color = Color::new(250, 128, 114);
    /// assert_eq!(color.green(), 128);
    /// ```
    pub const fn green(&self) -> u8 {
        self.g
    }

    /// Returns the value of the red component.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::feature::record_buf::Color;
    /// let color = Color::new(250, 128, 114);
    /// assert_eq!(color.blue(), 114);
    /// ```
    pub const fn blue(&self) -> u8 {
        self.b
    }
}

impl fmt::Display for Color {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}{}{}{}{}",
            self.r, DELIMITER, self.g, DELIMITER, self.b
        )
    }
}

/// An error returned when a raw BED record color fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
    /// A color component is invalid.
    InvalidComponent(num::ParseIntError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid => None,
            Self::InvalidComponent(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidComponent(e) => write!(f, "invalid component: {e}"),
        }
    }
}

impl FromStr for Color {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (a, b, c) = split_twice(s, DELIMITER).ok_or(ParseError::Invalid)?;

        let r = a.parse().map_err(ParseError::InvalidComponent)?;
        let g = b.parse().map_err(ParseError::InvalidComponent)?;
        let b = c.parse().map_err(ParseError::InvalidComponent)?;

        Ok(Color::new(r, g, b))
    }
}

fn split_twice(s: &str, delimiter: char) -> Option<(&str, &str, &str)> {
    let (a, s) = s.split_once(delimiter)?;
    let (b, c) = s.split_once(delimiter)?;
    Some((a, b, c))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Color::new(0, 0, 0).to_string(), "0,0,0");
        assert_eq!(Color::new(8, 16, 128).to_string(), "8,16,128");
        assert_eq!(Color::new(255, 255, 255).to_string(), "255,255,255");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("0,0,0".parse(), Ok(Color::new(0, 0, 0)));
        assert_eq!("8,16,128".parse(), Ok(Color::new(8, 16, 128)));
        assert_eq!("255,255,255".parse(), Ok(Color::new(255, 255, 255)));

        assert_eq!("0".parse::<Color>(), Err(ParseError::Invalid));
        assert_eq!("0,0".parse::<Color>(), Err(ParseError::Invalid));

        assert!(matches!(
            "0,0,0,0".parse::<Color>(),
            Err(ParseError::InvalidComponent(_))
        ));
        assert!(matches!(
            "red,0,0".parse::<Color>(),
            Err(ParseError::InvalidComponent(_))
        ));
        assert!(matches!(
            "0,green,0".parse::<Color>(),
            Err(ParseError::InvalidComponent(_))
        ));
        assert!(matches!(
            "0,0,blue".parse::<Color>(),
            Err(ParseError::InvalidComponent(_))
        ));
    }
}
