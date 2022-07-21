use std::str::FromStr;

use crate::header::record::value::map::tag;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Standard {
    Id,
    Type,
    Number,
    Values,
}

impl tag::Standard for Standard {}

impl FromStr for Standard {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "ID" => Ok(Self::Id),
            "Type" => Ok(Self::Type),
            "Number" => Ok(Self::Number),
            "Values" => Ok(Self::Values),
            _ => Err(()),
        }
    }
}
