mod context;
mod record;

use std::{error, fmt, hash::Hash, str};

use bstr::BString;
use indexmap::IndexMap;

pub(crate) use self::context::Context;
use self::record::parse_record;
use super::{
    Header, Programs, ReadGroups, Record, ReferenceSequences,
    record::value::{
        Map,
        map::{self, header::Version},
    },
};

/// An error returned when a raw SAM header fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A header record is not on the first line.
    UnexpectedHeader,
    /// The record is invalid.
    InvalidRecord(record::ParseError),
    /// A reference sequence name is duplicated.
    DuplicateReferenceSequenceName(BString),
    /// A read group ID is duplicated.
    DuplicateReadGroupId(BString),
    /// A program ID is duplicated.
    DuplicateProgramId(BString),
    /// A comment record is invalid.
    InvalidComment,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidRecord(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedHeader => write!(f, "unexpected header (HD) record"),
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

/// A SAM header parser.
#[derive(Default)]
pub struct Parser {
    ctx: Context,
    header: Option<Map<map::Header>>,
    reference_sequences: ReferenceSequences,
    read_groups: ReadGroups,
    programs: Programs,
    comments: Vec<BString>,
}

impl Parser {
    fn is_empty(&self) -> bool {
        self.header.is_none()
            && self.reference_sequences.is_empty()
            && self.read_groups.is_empty()
            && self.programs.as_ref().is_empty()
            && self.comments.is_empty()
    }

    /// Parses and adds a raw record to the header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let mut parser = sam::header::Parser::default();
    /// parser.parse_partial(b"@HD\tVN:1.6")?;
    /// # Ok::<_, sam::header::ParseError>(())
    /// ```
    pub fn parse_partial(&mut self, src: &[u8]) -> Result<(), ParseError> {
        if self.is_empty()
            && let Some(version) = extract_version(src)
        {
            self.ctx = Context::from(version);
        }

        let record = parse_record(src, &self.ctx).map_err(ParseError::InvalidRecord)?;

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
                self.programs.as_mut(),
                id,
                program,
                ParseError::DuplicateProgramId,
            )?,
            Record::Comment(comment) => self.comments.push(comment),
        }

        Ok(())
    }

    /// Builds the SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let parser = sam::header::Parser::default();
    /// let header = parser.finish();
    /// assert!(header.is_empty());
    /// # Ok::<_, sam::header::ParseError>(())
    /// ```
    pub fn finish(self) -> Header {
        Header {
            header: self.header,
            reference_sequences: self.reference_sequences,
            read_groups: self.read_groups,
            programs: self.programs,
            comments: self.comments,
        }
    }
}

fn extract_version(src: &[u8]) -> Option<Version> {
    use self::record::value::map::header::parse_version;

    const RECORD_PREFIX: &[u8] = b"@HD\t";
    const DELIMITER: u8 = b'\t';
    const FIELD_PREFIX: &[u8] = b"VN:";

    if let Some(raw_value) = src.strip_prefix(RECORD_PREFIX) {
        for raw_field in raw_value.split(|&b| b == DELIMITER) {
            if let Some(s) = raw_field.strip_prefix(FIELD_PREFIX) {
                return parse_version(s).ok();
            }
        }
    }

    None
}

fn try_insert<K, V, F, E>(map: &mut IndexMap<K, V>, key: K, value: V, f: F) -> Result<(), E>
where
    K: Hash + Eq + Clone,
    F: FnOnce(K) -> E,
{
    use indexmap::map::Entry;

    match map.entry(key) {
        Entry::Vacant(e) => {
            e.insert(value);
            Ok(())
        }
        Entry::Occupied(e) => Err(f(e.key().clone())),
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
/// assert!(header.programs().as_ref().is_empty());
/// assert!(header.comments().is_empty());
/// # Ok::<(), sam::header::ParseError>(())
/// ```
pub(super) fn parse(s: &str) -> Result<Header, ParseError> {
    let mut parser = Parser::default();

    for line in s.lines() {
        parser.parse_partial(line.as_bytes())?;
    }

    Ok(parser.finish())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZero;

        use crate::header::record::value::map::{
            self, Map, Program, ReadGroup, ReferenceSequence,
            header::{self, Version},
            program,
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
                    .insert(header::tag::SORT_ORDER, "coordinate")
                    .build()?,
            )
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
            )
            .add_reference_sequence(
                "sq1",
                Map::<ReferenceSequence>::new(const { NonZero::new(13).unwrap() }),
            )
            .add_read_group("rg0", Map::<ReadGroup>::default())
            .add_program(
                "pg0",
                Map::<Program>::builder()
                    .insert(program::tag::NAME, "noodles")
                    .build()?,
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
        assert!(header.programs().as_ref().is_empty());
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
    fn test_parse_with_duplicate_reference_sequence_names() {
        let s = "\
@SQ\tSN:sq0\tLN:8
@SQ\tSN:sq0\tLN:8
";

        assert_eq!(
            parse(s),
            Err(ParseError::DuplicateReferenceSequenceName(BString::from(
                "sq0"
            )))
        );
    }

    #[test]
    fn test_parse_with_duplicate_read_group_ids() {
        let s = "\
@RG\tID:rg0
@RG\tID:rg0
";

        assert_eq!(
            parse(s),
            Err(ParseError::DuplicateReadGroupId(BString::from("rg0")))
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
            Err(ParseError::DuplicateProgramId(BString::from("pg0")))
        );
    }

    #[test]
    fn test_extract_version() {
        assert_eq!(extract_version(b"@HD\tVN:1.6"), Some(Version::new(1, 6)));
        assert_eq!(
            extract_version(b"@HD\tSO:coordinate\tVN:1.6"),
            Some(Version::new(1, 6))
        );
        assert!(extract_version(b"@HD\tVN:NA").is_none());
        assert!(extract_version(b"@SQ\tSN:sq0\tLN:8\tVN:1.6").is_none());
        assert!(extract_version(b"@CO\tVN:1.6").is_none());
    }
}
