use std::fmt::{self, Write};

use indexmap::IndexMap;

/// A VCF header record value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Value {
    /// A string.
    String(String),
    /// A structure.
    Struct(String, IndexMap<String, String>),
}

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::String(value) => f.write_str(value),
            Self::Struct(id, fields) => {
                f.write_char('<')?;

                write!(f, "ID={}", id)?;

                for (key, value) in fields {
                    write!(f, ",{}={}", key, value)?;
                }

                f.write_char('>')?;

                Ok(())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(
            Value::String(String::from("VCFv4.3")).to_string(),
            "VCFv4.3"
        );

        assert_eq!(
            Value::Struct(String::from("sq0"), IndexMap::new()).to_string(),
            "<ID=sq0>"
        );

        assert_eq!(
            Value::Struct(
                String::from("sq0"),
                [(String::from("length"), String::from("13"))]
                    .into_iter()
                    .collect(),
            )
            .to_string(),
            "<ID=sq0,length=13>"
        );
    }
}
