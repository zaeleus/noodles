mod context;
mod record;

use std::{error, fmt, hash::Hash};

use indexmap::IndexMap;

pub(crate) use self::context::Context;
use super::{
    record::value::{
        map::{self, reference_sequence},
        Map,
    },
    Header, Programs, ReadGroups, Record, ReferenceSequences,
};

/// An error returned when a raw SAM header fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A header record is not on the first line.
    UnexpectedHeader,
    /// The header record is invalid.
    InvalidHeaderRecord(super::record::ParseError),
    /// The record is invalid.
    InvalidRecord(record::ParseError),
    /// A reference sequence name is duplicated.
    DuplicateReferenceSequenceName(reference_sequence::Name),
    /// A read group ID is duplicated.
    DuplicateReadGroupId(String),
    /// A program ID is duplicated.
    DuplicateProgramId(String),
    /// A comment record is invalid.
    InvalidComment,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidHeaderRecord(e) => Some(e),
            Self::InvalidRecord(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedHeader => f.write_str("unexpected @HD"),
            Self::InvalidHeaderRecord(_) => f.write_str("invalid header record"),
            Self::InvalidRecord(_) => f.write_str("invalid record"),
            Self::DuplicateReferenceSequenceName(name) => {
                write!(f, "duplicate reference sequence name: {name}")
            }
            Self::DuplicateReadGroupId(id) => write!(f, "duplicate read group ID: {id}"),
            Self::DuplicateProgramId(id) => write!(f, "duplicate program ID: {id}"),
            Self::InvalidComment => f.write_str("invalid comment record"),
        }
    }
}

#[derive(Default)]
pub(crate) struct Parser {
    ctx: Context,
    header: Option<Map<map::Header>>,
    reference_sequences: ReferenceSequences,
    read_groups: ReadGroups,
    programs: Programs,
    comments: Vec<String>,
}

impl Parser {
    fn is_empty(&self) -> bool {
        self.header.is_none()
            && self.reference_sequences.is_empty()
            && self.read_groups.is_empty()
            && self.programs.is_empty()
            && self.comments.is_empty()
    }

    pub(crate) fn parse_partial(&mut self, s: &str) -> Result<(), ParseError> {
        if self.is_empty() {
            if let Some(result) = super::record::extract_version(s) {
                let version = result.map_err(ParseError::InvalidHeaderRecord)?;
                self.ctx = Context::from(version);
            }
        }

        let record =
            record::parse_record(s.as_bytes(), &self.ctx).map_err(ParseError::InvalidRecord)?;

        match record {
            Record::Header(header) => {
                if self.is_empty() {
                    self.header = Some(header);
                } else {
                    return Err(ParseError::UnexpectedHeader);
                }
            }
            Record::ReferenceSequence(name, reference_sequence) => try_insert(
                &mut self.reference_sequences,
                name,
                reference_sequence,
                ParseError::DuplicateReferenceSequenceName,
            )?,
            Record::ReadGroup(id, read_group) => try_insert(
                &mut self.read_groups,
                id,
                read_group,
                ParseError::DuplicateReadGroupId,
            )?,
            Record::Program(id, program) => try_insert(
                &mut self.programs,
                id,
                program,
                ParseError::DuplicateProgramId,
            )?,
            Record::Comment(comment) => self.comments.push(comment),
        }

        Ok(())
    }

    pub(crate) fn finish(self) -> Header {
        Header {
            header: self.header,
            reference_sequences: self.reference_sequences,
            read_groups: self.read_groups,
            programs: self.programs,
            comments: self.comments,
        }
    }
}

fn try_insert<K, V, F, E>(map: &mut IndexMap<K, V>, key: K, value: V, f: F) -> Result<(), E>
where
    K: Hash + Eq,
    F: FnOnce(K) -> E,
{
    use indexmap::map::Entry;

    match map.entry(key) {
        Entry::Vacant(e) => {
            e.insert(value);
            Ok(())
        }
        Entry::Occupied(e) => {
            let (k, _) = e.remove_entry();
            Err(f(k))
        }
    }
}

/// Parses a raw SAM header.
///
/// # Examples
///
/// ```
/// use noodles_sam as sam;
///
/// let s = "\
/// @HD\tVN:1.6\tSO:coordinate
/// @SQ\tSN:sq0\tLN:8
/// @SQ\tSN:sq1\tLN:13
/// ";
///
/// let header: sam::Header = s.parse()?;
///
/// assert!(header.header().is_some());
/// assert_eq!(header.reference_sequences().len(), 2);
/// assert!(header.read_groups().is_empty());
/// assert!(header.programs().is_empty());
/// assert!(header.comments().is_empty());
/// # Ok::<(), sam::header::ParseError>(())
/// ```
pub(super) fn parse(s: &str) -> Result<Header, ParseError> {
    let mut parser = Parser::default();

    for line in s.lines() {
        parser.parse_partial(line)?;
    }

    Ok(parser.finish())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZeroUsize;

        use crate::header::record::value::map::{
            self,
            header::{SortOrder, Version},
            Map, Program, ReadGroup, ReferenceSequence,
        };

        let s = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:sq0\tLN:8
@SQ\tSN:sq1\tLN:13
@RG\tID:rg0
@PG\tID:pg0\tPN:noodles
@CO\tndls
";

        let actual = parse(s)?;

        let expected = Header::builder()
            .set_header(
                Map::<map::Header>::builder()
                    .set_version(Version::new(1, 6))
                    .set_sort_order(SortOrder::Coordinate)
                    .build()?,
            )
            .add_reference_sequence(
                "sq0".parse()?,
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8)?),
            )
            .add_reference_sequence(
                "sq1".parse()?,
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?),
            )
            .add_read_group("rg0", Map::<ReadGroup>::default())
            .add_program(
                "pg0",
                Map::<Program>::builder().set_name("noodles").build()?,
            )
            .add_comment("ndls")
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_parse_with_empty_input() -> Result<(), ParseError> {
        let header = parse("")?;

        assert!(header.header().is_none());
        assert!(header.reference_sequences().is_empty());
        assert!(header.read_groups().is_empty());
        assert!(header.programs().is_empty());
        assert!(header.comments().is_empty());

        Ok(())
    }

    #[test]
    fn test_parse_without_hd() -> Result<(), ParseError> {
        let header = parse("@SQ\tSN:sq0\tLN:8\n")?;
        assert!(header.header().is_none());
        assert_eq!(header.reference_sequences().len(), 1);
        Ok(())
    }

    #[test]
    fn test_parse_with_multiple_hd() {
        let s = "\
@HD\tVN:1.6\tSO:coordinate
@HD\tVN:1.6\tSO:coordinate
";

        assert_eq!(parse(s), Err(ParseError::UnexpectedHeader));
    }

    #[test]
    fn test_parse_with_duplicate_reference_sequence_names(
    ) -> Result<(), reference_sequence::name::ParseError> {
        let s = "\
@SQ\tSN:sq0\tLN:8
@SQ\tSN:sq0\tLN:8
";

        assert_eq!(
            parse(s),
            Err(ParseError::DuplicateReferenceSequenceName("sq0".parse()?))
        );

        Ok(())
    }

    #[test]
    fn test_parse_with_duplicate_read_group_ids() {
        let s = "\
@RG\tID:rg0
@RG\tID:rg0
";

        assert_eq!(
            parse(s),
            Err(ParseError::DuplicateReadGroupId(String::from("rg0")))
        );
    }

    #[test]
    fn test_parse_with_duplicate_program_ids() {
        let s = "\
@PG\tID:pg0
@PG\tID:pg0
";
        assert_eq!(
            parse(s),
            Err(ParseError::DuplicateProgramId(String::from("pg0")))
        );
    }
}
