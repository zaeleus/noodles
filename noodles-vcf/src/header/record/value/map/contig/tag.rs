use std::str::FromStr;

use crate::header::record::value::map::tag;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Standard {
    Id,
    Length,
    Idx,
}

impl tag::Standard for Standard {}

impl FromStr for Standard {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "ID" => Ok(Self::Id),
            "length" => Ok(Self::Length),
            "IDX" => Ok(Self::Idx),
            _ => Err(()),
        }
    }
}
