use std::{borrow::Cow, io, iter};

use crate::io::reader::record_buf::value::percent_decode;

/// Variant record sample array value values.
#[allow(clippy::len_without_is_empty)]
pub trait Values<'a, N> {
    /// Returns the number of values.
    fn len(&self) -> usize;

    /// Returns an iterator over values.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<N>>> + '_>;
}

const DELIMITER: char = ',';
const MISSING: &str = ".";

impl<'a> Values<'a, i32> for &'a str {
    fn len(&self) -> usize {
        count(self)
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<i32>>> + '_> {
        if self.is_empty() {
            Box::new(iter::empty())
        } else {
            Box::new(self.split(DELIMITER).map(|s| {
                match s {
                    MISSING => Ok(None),
                    _ => s
                        .parse()
                        .map(Some)
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
                }
            }))
        }
    }
}

impl<'a> Values<'a, f32> for &'a str {
    fn len(&self) -> usize {
        count(self)
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<f32>>> + '_> {
        if self.is_empty() {
            Box::new(iter::empty())
        } else {
            Box::new(self.split(DELIMITER).map(|s| {
                match s {
                    MISSING => Ok(None),
                    _ => s
                        .parse()
                        .map(Some)
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
                }
            }))
        }
    }
}

impl<'a> Values<'a, char> for &'a str {
    fn len(&self) -> usize {
        count(self)
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<char>>> + '_> {
        const MISSING: &str = ".";

        if self.is_empty() {
            Box::new(iter::empty())
        } else {
            Box::new(self.split(DELIMITER).map(|s| {
                match s {
                    MISSING => Ok(None),
                    _ => percent_decode(s)
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                        .and_then(|t| {
                            let mut chars = t.chars();

                            if let Some(c) = chars.next()
                                && chars.next().is_none()
                            {
                                return Ok(Some(c));
                            }

                            Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "invalid character",
                            ))
                        }),
                }
            }))
        }
    }
}

impl<'a> Values<'a, Cow<'a, str>> for &'a str {
    fn len(&self) -> usize {
        count(self)
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<Cow<'a, str>>>> + '_> {
        if self.is_empty() {
            Box::new(iter::empty())
        } else {
            Box::new(self.split(DELIMITER).map(|s| {
                match s {
                    MISSING => Ok(None),
                    _ => percent_decode(s)
                        .map(Some)
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
                }
            }))
        }
    }
}

fn count(s: &str) -> usize {
    const DELIMITER: u8 = b',';

    if s.is_empty() {
        0
    } else {
        let n = s.as_bytes().iter().filter(|&&b| b == DELIMITER).count();
        n + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integer_values() -> io::Result<()> {
        let src = "";
        let values: Box<dyn Values<'_, i32>> = Box::new(src);
        assert_eq!(values.len(), 0);

        let mut iter = values.iter();
        assert!(iter.next().transpose()?.is_none());

        let src = "8,13,.";
        let values: Box<dyn Values<'_, i32>> = Box::new(src);
        assert_eq!(values.len(), 3);

        let mut iter = values.iter();
        assert_eq!(iter.next().transpose()?, Some(Some(8)));
        assert_eq!(iter.next().transpose()?, Some(Some(13)));
        assert_eq!(iter.next().transpose()?, Some(None));
        assert!(iter.next().transpose()?.is_none());

        Ok(())
    }

    #[test]
    fn test_float_values() -> io::Result<()> {
        let src = "";
        let values: Box<dyn Values<'_, f32>> = Box::new(src);
        assert_eq!(values.len(), 0);

        let mut iter = values.iter();
        assert!(iter.next().transpose()?.is_none());

        let src = "8.0,13.0,.";
        let values: Box<dyn Values<'_, f32>> = Box::new(src);
        assert_eq!(values.len(), 3);

        let mut iter = values.iter();
        assert_eq!(iter.next().transpose()?, Some(Some(8.0)));
        assert_eq!(iter.next().transpose()?, Some(Some(13.0)));
        assert_eq!(iter.next().transpose()?, Some(None));
        assert!(iter.next().transpose()?.is_none());

        Ok(())
    }

    #[test]
    fn test_character_values() -> io::Result<()> {
        let src = "";
        let values: Box<dyn Values<'_, char>> = Box::new(src);
        assert_eq!(values.len(), 0);

        let mut iter = values.iter();
        assert!(iter.next().transpose()?.is_none());

        let src = "a,%3A,.";
        let values: Box<dyn Values<'_, char>> = Box::new(src);
        assert_eq!(values.len(), 3);

        let mut iter = values.iter();
        assert_eq!(iter.next().transpose()?, Some(Some('a')));
        assert_eq!(iter.next().transpose()?, Some(Some(':')));
        assert_eq!(iter.next().transpose()?, Some(None));
        assert!(iter.next().transpose()?.is_none());

        Ok(())
    }

    #[test]
    fn test_string_values() -> io::Result<()> {
        let src = "";
        let values: Box<dyn Values<'_, Cow<'_, str>>> = Box::new(src);
        assert_eq!(values.len(), 0);

        let mut iter = values.iter();
        assert!(iter.next().transpose()?.is_none());

        let src = "a,b%3Ac,.";
        let values: Box<dyn Values<'_, Cow<'_, str>>> = Box::new(src);
        assert_eq!(values.len(), 3);

        let mut iter = values.iter();
        assert_eq!(iter.next().transpose()?, Some(Some(Cow::from("a"))));
        assert_eq!(iter.next().transpose()?, Some(Some(Cow::from("b:c"))));
        assert_eq!(iter.next().transpose()?, Some(None));
        assert!(iter.next().transpose()?.is_none());

        Ok(())
    }
}
