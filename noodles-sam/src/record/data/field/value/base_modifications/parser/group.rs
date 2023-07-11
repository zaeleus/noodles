use std::{error, fmt};

use crate::record::{
    data::field::value::base_modifications::{group::UnmodifiedBase, Group},
    Sequence,
};

mod modifications;
mod status;
mod strand;
mod unmodified_base;

use self::{
    modifications::parse_modifications, status::parse_status, strand::parse_strand,
    unmodified_base::parse_unmodified_base,
};

/// An error returned when a base modifications group fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The unmodified base is invalid.
    InvalidUnmodifiedBase(unmodified_base::ParseError),
    /// The strand is invalid.
    InvalidStrand(strand::ParseError),
    /// A modification is invalid.
    InvalidModifications(modifications::ParseError),
    /// The status is invalid.
    InvalidStatus,
    /// A skip count is invalid.
    InvalidSkipCount(lexical_core::Error),
    /// The terminator is invalid.
    InvalidTerminator,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidUnmodifiedBase(e) => Some(e),
            Self::InvalidStrand(e) => Some(e),
            Self::InvalidModifications(e) => Some(e),
            Self::InvalidSkipCount(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidUnmodifiedBase(_) => write!(f, "invalid unmodified base"),
            Self::InvalidStrand(_) => write!(f, "invalid strand"),
            Self::InvalidModifications(_) => write!(f, "invalid modifications"),
            Self::InvalidStatus => write!(f, "invalid status"),
            Self::InvalidSkipCount(_) => write!(f, "invalid skip count"),
            Self::InvalidTerminator => write!(f, "invalid terminator"),
        }
    }
}

pub(super) fn parse_group(src: &mut &[u8], sequence: &Sequence) -> Result<Group, ParseError> {
    let unmodified_base = parse_unmodified_base(src).map_err(ParseError::InvalidUnmodifiedBase)?;
    let strand = parse_strand(src).map_err(ParseError::InvalidStrand)?;
    let modifications = parse_modifications(src).map_err(ParseError::InvalidModifications)?;
    let status = parse_status(src);
    let skip_counts = parse_skip_counts(src)?;
    consume_terminator(src)?;

    let positions = decode_positions(&skip_counts, sequence, unmodified_base)?;

    Ok(Group::new(
        unmodified_base,
        strand,
        modifications,
        status,
        positions,
    ))
}

fn parse_skip_counts(src: &mut &[u8]) -> Result<Vec<usize>, ParseError> {
    const DELIMITER: u8 = b',';

    let mut counts = Vec::new();

    loop {
        if let Some((b, rest)) = src.split_first() {
            if *b == DELIMITER {
                *src = rest;
            } else {
                break;
            }
        } else {
            return Err(ParseError::UnexpectedEof);
        }

        let (n, i) = lexical_core::parse_partial(src).map_err(ParseError::InvalidSkipCount)?;
        *src = &src[i..];
        counts.push(n);
    }

    Ok(counts)
}

fn consume_terminator(src: &mut &[u8]) -> Result<(), ParseError> {
    const TERMINATOR: u8 = b';';

    if let Some((b, rest)) = src.split_first() {
        if *b == TERMINATOR {
            *src = rest;
            Ok(())
        } else {
            Err(ParseError::InvalidTerminator)
        }
    } else {
        Err(ParseError::UnexpectedEof)
    }
}

fn decode_positions(
    skip_counts: &[usize],
    sequence: &Sequence,
    unmodified_base: UnmodifiedBase,
) -> Result<Vec<usize>, ParseError> {
    use crate::record::sequence::Base;

    let mut positions = Vec::with_capacity(skip_counts.len());
    let unmodified_base = Base::from(unmodified_base);

    let mut iter = sequence
        .as_ref()
        .iter()
        .enumerate()
        .filter(|(_, &base)| base == unmodified_base)
        .map(|(i, _)| i);

    for &count in skip_counts {
        let i = iter.nth(count).unwrap();
        positions.push(i);
    }

    Ok(positions)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_group() -> Result<(), crate::record::sequence::ParseError> {
        use crate::record::data::field::value::base_modifications::group::{
            modification, Status, Strand, UnmodifiedBase,
        };

        let sequence = "CACCCGATGACCGGCT".parse()?;

        let mut src = &b"C+m,1,3,0;"[..];
        let actual = parse_group(&mut src, &sequence);
        let expected = Group::new(
            UnmodifiedBase::C,
            Strand::Forward,
            vec![modification::FIVE_METHYLCYTOSINE],
            None,
            vec![2, 11, 14],
        );
        assert_eq!(actual, Ok(expected));

        let mut src = &b"C+m.,1,3,0;"[..];
        let actual = parse_group(&mut src, &sequence);
        let expected = Group::new(
            UnmodifiedBase::C,
            Strand::Forward,
            vec![modification::FIVE_METHYLCYTOSINE],
            Some(Status::Implicit),
            vec![2, 11, 14],
        );
        assert_eq!(actual, Ok(expected));

        let mut src = &b""[..];
        assert!(matches!(
            parse_group(&mut src, &sequence),
            Err(ParseError::InvalidUnmodifiedBase(_))
        ));

        let mut src = &b"C"[..];
        assert!(matches!(
            parse_group(&mut src, &sequence),
            Err(ParseError::InvalidStrand(_))
        ));

        let mut src = &b"C+"[..];
        assert!(matches!(
            parse_group(&mut src, &sequence),
            Err(ParseError::InvalidModifications(_))
        ));

        let mut src = &b"C+m,"[..];
        assert!(matches!(
            parse_group(&mut src, &sequence),
            Err(ParseError::InvalidSkipCount(_))
        ));

        Ok(())
    }

    #[test]
    fn test_consume_terminator() {
        let mut src = &b";"[..];
        assert!(consume_terminator(&mut src).is_ok());

        let mut src = &b"n"[..];
        assert_eq!(
            consume_terminator(&mut src),
            Err(ParseError::InvalidTerminator)
        );

        let mut src = &b""[..];
        assert_eq!(consume_terminator(&mut src), Err(ParseError::UnexpectedEof));
    }
}
