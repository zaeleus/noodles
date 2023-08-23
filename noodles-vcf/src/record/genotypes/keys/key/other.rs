use std::{fmt, str::FromStr};

use super::ParseError;
use crate::header::FileFormat;

/// A non-reserved VCF header format key.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Other(pub(super) String);

impl AsRef<str> for Other {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for Other {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.as_ref().fmt(f)
    }
}

impl FromStr for Other {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if is_valid_name(s) {
            Ok(Self(s.into()))
        } else {
            Err(ParseError::Invalid)
        }
    }
}

impl TryFrom<(FileFormat, &str)> for Other {
    type Error = ParseError;

    fn try_from((file_format, s): (FileFormat, &str)) -> Result<Self, Self::Error> {
        const VCF_4_3: FileFormat = FileFormat::new(4, 3);

        if file_format < VCF_4_3 {
            Ok(Self(s.into()))
        } else {
            s.parse()
        }
    }
}

// § 1.6.2 Genotype fields
fn is_valid_name_char(c: char) -> bool {
    matches!(c, '0'..='9' | 'A'..='Z' | 'a'..='z' | '_' | '.')
}

fn is_valid_name(s: &str) -> bool {
    let mut chars = s.chars();

    if let Some(c) = chars.next() {
        if !matches!(c, 'A'..='Z' | 'a'..='z' | '_') {
            return false;
        }
    }

    chars.all(is_valid_name_char)
}
