use std::fmt;

const DELIMITER: char = ',';
const MISSING: char = '.';

/// A VCF record info field array value.
#[derive(Clone, Debug, PartialEq)]
pub enum Array {
    /// An array of 32-bit integers.
    Integer(Vec<Option<i32>>),
    /// An array of single-precision floating-points.
    Float(Vec<Option<f32>>),
    /// An array of characters.
    Character(Vec<Option<char>>),
    /// An array of strings.
    String(Vec<Option<String>>),
}

impl fmt::Display for Array {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Integer(values) => {
                for (i, value) in values.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{DELIMITER}")?;
                    }

                    if let Some(v) = value {
                        write!(f, "{v}")?;
                    } else {
                        write!(f, "{MISSING}")?;
                    }
                }

                Ok(())
            }
            Self::Float(values) => {
                for (i, value) in values.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{DELIMITER}")?;
                    }

                    if let Some(v) = value {
                        write!(f, "{v}")?;
                    } else {
                        write!(f, "{MISSING}")?;
                    }
                }

                Ok(())
            }
            Self::Character(values) => {
                for (i, value) in values.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{DELIMITER}")?;
                    }

                    if let Some(v) = value {
                        write!(f, "{v}")?;
                    } else {
                        write!(f, "{MISSING}")?;
                    }
                }

                Ok(())
            }
            Self::String(values) => {
                for (i, value) in values.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{DELIMITER}")?;
                    }

                    if let Some(v) = value {
                        write!(f, "{v}")?;
                    } else {
                        write!(f, "{MISSING}")?;
                    }
                }

                Ok(())
            }
        }
    }
}
