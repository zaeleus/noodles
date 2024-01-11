mod group;

use std::{error, fmt};

use self::group::parse_group;
use super::BaseModifications;
use crate::alignment::record_buf::Sequence;

/// An error returned when base modifications fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidGroup(group::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidGroup(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidGroup(_) => write!(f, "invalid group"),
        }
    }
}

pub(super) fn parse(
    s: &str,
    is_reverse_complemented: bool,
    sequence: &Sequence,
) -> Result<BaseModifications, ParseError> {
    let mut groups = Vec::new();
    let mut src = s.as_bytes();

    while !src.is_empty() {
        let group = parse_group(&mut src, is_reverse_complemented, sequence)
            .map_err(ParseError::InvalidGroup)?;

        groups.push(group);
    }

    Ok(BaseModifications(groups))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() {
        use crate::record::data::field::value::base_modifications::{
            group::{modification, Strand, UnmodifiedBase},
            Group,
        };

        let is_reverse_complemented = false;
        let sequence = Sequence::from(b"CACCCGATGACCGGCT".to_vec());
        let actual = parse("C+m,1,3,0;G-o,2;", is_reverse_complemented, &sequence);

        let expected = BaseModifications(vec![
            Group::new(
                UnmodifiedBase::C,
                Strand::Forward,
                vec![modification::FIVE_METHYLCYTOSINE],
                None,
                vec![2, 11, 14],
            ),
            Group::new(
                UnmodifiedBase::G,
                Strand::Reverse,
                vec![modification::EIGHT_OXOGUANINE],
                None,
                vec![12],
            ),
        ]);

        assert_eq!(actual, Ok(expected));
    }
}
