//! SAM header and records.
//!
//! A SAM header defines 5 record types:
//!
//!   1. header (`@HD`),
//!   2. reference sequence (`@SQ`),
//!   3. read group (`@RG`),
//!   4. program (`@PG`), and
//!   5. comment (`@CO`).
//!
//! Each record is effectively a map. It defines key-value pairs associated with that record type.
//!
//! All records are optional, which means an empty header is considered a valid SAM header.
//!
//! If there is a header record, it must appear on the first line.
//!
//! Reference sequence, read group, program, and comment records are lists of records of the same
//! type. Reference sequences must be ordered; whereas read groups, programs, and comments can be
//! unordered. (`sam::Header` defines them to be ordered.)
//!
//! # Examples
//!
//! ## Parse a SAM header
//!
//! ```
//! use noodles_sam as sam;
//!
//! let s = "\
//! @HD\tVN:1.6\tSO:coordinate
//! @SQ\tSN:sq0\tLN:8
//! @SQ\tSN:sq1\tLN:13
//! ";
//!
//! let header: sam::Header = s.parse()?;
//!
//! assert!(header.header().is_some());
//! assert_eq!(header.reference_sequences().len(), 2);
//! assert!(header.read_groups().is_empty());
//! assert!(header.programs().is_empty());
//! assert!(header.comments().is_empty());
//! # Ok::<(), sam::header::ParseError>(())
//! ```
//!
//! ## Create a SAM header
//!
//! ```
//! use noodles_sam::{self as sam, header::{self, ReferenceSequence}};
//!
//! let header = sam::Header::builder()
//!     .set_header(header::header::Header::default())
//!     .add_reference_sequence(ReferenceSequence::new(String::from("sq0"), 8))
//!     .add_reference_sequence(ReferenceSequence::new(String::from("sq1"), 13))
//!     .build();
//!
//! assert!(header.header().is_some());
//! assert_eq!(header.reference_sequences().len(), 2);
//! assert!(header.read_groups().is_empty());
//! assert!(header.programs().is_empty());
//! assert!(header.comments().is_empty());
//! ```

mod builder;
#[allow(clippy::module_inception)]
pub mod header;
pub mod program;
pub mod read_group;
pub mod record;
pub mod reference_sequence;

use std::{convert::TryFrom, error, fmt, str::FromStr};

use indexmap::IndexMap;

pub use self::{
    builder::Builder, program::Program, read_group::ReadGroup,
    reference_sequence::ReferenceSequence,
};

pub use self::record::Record;

/// A reference seqeuence dictionary.
pub type ReferenceSequences = IndexMap<String, ReferenceSequence>;

/// An ordered map of read groups.
pub type ReadGroups = IndexMap<String, ReadGroup>;

/// An ordered map of programs.
pub type Programs = IndexMap<String, Program>;

/// A SAM header.
///
/// Records are grouped by their types: header, reference seqeuence, read group, program, and
/// comment.
#[derive(Debug, Default)]
pub struct Header {
    header: Option<header::Header>,
    reference_sequences: ReferenceSequences,
    read_groups: ReadGroups,
    programs: Programs,
    comments: Vec<String>,
}

impl Header {
    /// Returns a builder to create a SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::{self, ReferenceSequence}};
    ///
    /// let header = sam::Header::builder()
    ///     .set_header(header::header::Header::default())
    ///     .add_reference_sequence(ReferenceSequence::new(String::from("sq0"), 13))
    ///     .build();
    ///
    /// assert!(header.header().is_some());
    /// assert_eq!(header.reference_sequences().len(), 1);
    /// assert!(header.read_groups().is_empty());
    /// assert!(header.programs().is_empty());
    /// assert!(header.comments().is_empty());
    /// ```
    pub fn builder() -> Builder {
        Builder::new()
    }

    /// Returns the SAM header header if it is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header};
    ///
    /// let header = sam::Header::default();
    /// assert!(header.header().is_none());
    ///
    /// let header = sam::Header::builder()
    ///     .set_header(header::header::Header::default())
    ///     .build();
    ///
    /// assert!(header.header().is_some());
    /// ```
    pub fn header(&self) -> Option<&header::Header> {
        self.header.as_ref()
    }

    /// Returns a mutable reference to the SAM header header if it is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header};
    ///
    /// let mut header = sam::Header::builder()
    ///     .set_header(header::header::Header::new(String::from("1.6")))
    ///     .build();
    /// assert_eq!(header.header().map(|h| h.version()), Some("1.6"));
    ///
    /// header.header_mut().as_mut().map(|h| {
    ///     *h.version_mut() = String::from("1.5");
    /// });
    /// assert_eq!(header.header().map(|h| h.version()), Some("1.5"));
    ///
    /// *header.header_mut() = None;
    /// assert!(header.header().is_none());
    /// ```
    pub fn header_mut(&mut self) -> &mut Option<header::Header> {
        &mut self.header
    }

    /// Returns the SAM header reference sequences.
    ///
    /// This is also called the reference sequence dictionary.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::ReferenceSequence};
    ///
    /// let header = sam::Header::builder()
    ///     .add_reference_sequence(ReferenceSequence::new(String::from("sq0"), 13))
    ///     .build();
    ///
    /// let reference_sequences = header.reference_sequences();
    /// assert_eq!(reference_sequences.len(), 1);
    /// assert!(reference_sequences.contains_key("sq0"));
    /// ```
    pub fn reference_sequences(&self) -> &ReferenceSequences {
        &self.reference_sequences
    }

    /// Returns a mutable reference to the SAM header reference sequences.
    ///
    /// This is also called the reference sequence dictionary.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::ReferenceSequence};
    ///
    /// let mut header = sam::Header::default();
    /// assert!(header.reference_sequences().is_empty());
    ///
    /// header.reference_sequences_mut().insert(
    ///     String::from("sq0"),
    ///     ReferenceSequence::new(String::from("sq0"), 13)
    /// );
    ///
    /// let reference_sequences = header.reference_sequences();
    /// assert_eq!(reference_sequences.len(), 1);
    /// assert!(reference_sequences.contains_key("sq0"));
    /// ```
    pub fn reference_sequences_mut(&mut self) -> &mut ReferenceSequences {
        &mut self.reference_sequences
    }

    /// Returns the SAM header read groups.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::ReadGroup};
    ///
    /// let header = sam::Header::builder()
    ///     .add_read_group(ReadGroup::new(String::from("rg0")))
    ///     .build();
    ///
    /// let read_groups = header.read_groups();
    /// assert_eq!(read_groups.len(), 1);
    /// assert!(read_groups.contains_key("rg0"));
    /// ```
    pub fn read_groups(&self) -> &ReadGroups {
        &self.read_groups
    }

    /// Returns a mutable reference to the SAM header read groups.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::ReadGroup};
    ///
    /// let mut header = sam::Header::default();
    /// assert!(header.read_groups().is_empty());
    ///
    /// header.read_groups_mut().insert(
    ///     String::from("rg0"),
    ///     ReadGroup::new(String::from("rg0")),
    /// );
    ///
    /// let read_groups = header.read_groups();
    /// assert_eq!(read_groups.len(), 1);
    /// assert!(read_groups.contains_key("rg0"));
    /// ```
    pub fn read_groups_mut(&mut self) -> &mut ReadGroups {
        &mut self.read_groups
    }

    /// Returns the SAM header programs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::Program};
    ///
    /// let header = sam::Header::builder()
    ///     .add_program(Program::new(String::from("noodles-sam")))
    ///     .build();
    ///
    /// let programs = header.programs();
    /// assert_eq!(programs.len(), 1);
    /// assert!(programs.contains_key("noodles-sam"));
    /// ```
    pub fn programs(&self) -> &Programs {
        &self.programs
    }

    /// Returns a mutable reference to the SAM header programs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::Program};
    ///
    /// let mut header = sam::Header::default();
    ///
    /// header.programs_mut().insert(
    ///     String::from("noodles-sam"),
    ///     Program::new(String::from("noodles-sam")),
    /// );
    ///
    /// let programs = header.programs();
    /// assert_eq!(programs.len(), 1);
    /// assert!(programs.contains_key("noodles-sam"));
    /// ```
    pub fn programs_mut(&mut self) -> &mut Programs {
        &mut self.programs
    }

    /// Returns the SAM header comments.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let header = sam::Header::builder().add_comment("noodles-sam").build();
    /// let comments = header.comments();
    /// assert_eq!(comments.len(), 1);
    /// assert_eq!(&comments[0], "noodles-sam");
    /// ```
    pub fn comments(&self) -> &[String] {
        &self.comments
    }

    /// Returns a mutable reference to the SAM header comments.
    ///
    /// To simply append a comment record, consider using [`Self::add_comment`] instead.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let mut header = sam::Header::default();
    /// header.comments_mut().push(String::from("noodles-sam"));
    ///
    /// let comments = header.comments();
    /// assert_eq!(comments.len(), 1);
    /// assert_eq!(&comments[0], "noodles-sam");
    /// ```
    pub fn comments_mut(&mut self) -> &mut Vec<String> {
        &mut self.comments
    }

    /// Adds a comment.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let mut header = sam::Header::default();
    /// header.add_comment("noodles-sam");
    ///
    /// let comments = header.comments();
    /// assert_eq!(comments.len(), 1);
    /// assert_eq!(&comments[0], "noodles-sam");
    /// ```
    pub fn add_comment<S>(&mut self, comment: S)
    where
        S: Into<String>,
    {
        self.comments.push(comment.into());
    }

    /// Returns whether there are no records in this SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let header = sam::Header::default();
    /// assert!(header.is_empty());
    ///
    /// let header = sam::Header::builder().add_comment("noodles-sam").build();
    /// assert!(!header.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.header.is_none()
            && self.reference_sequences.is_empty()
            && self.read_groups.is_empty()
            && self.programs.is_empty()
            && self.comments.is_empty()
    }
}

impl fmt::Display for Header {
    /// Formats the SAM header as a raw SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::{self, ReferenceSequence}};
    ///
    /// let header = sam::Header::builder()
    ///     .set_header(header::header::Header::new(String::from("1.6")))
    ///     .add_reference_sequence(ReferenceSequence::new(String::from("sq0"), 8))
    ///     .add_reference_sequence(ReferenceSequence::new(String::from("sq1"), 13))
    ///     .build();
    ///
    /// let expected = "\
    /// @HD\tVN:1.6
    /// @SQ\tSN:sq0\tLN:8
    /// @SQ\tSN:sq1\tLN:13
    /// ";
    ///
    /// assert_eq!(header.to_string(), expected);
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(header) = self.header() {
            writeln!(f, "{}", header)?;
        }

        for reference_sequence in self.reference_sequences.values() {
            writeln!(f, "{}", reference_sequence)?;
        }

        for read_group in self.read_groups.values() {
            writeln!(f, "{}", read_group)?;
        }

        for program in self.programs.values() {
            writeln!(f, "{}", program)?;
        }

        for comment in &self.comments {
            writeln!(f, "{}\t{}", record::Kind::Comment, comment)?;
        }

        Ok(())
    }
}

/// An error returned when a raw SAM header fails to parse.
#[derive(Debug)]
pub enum ParseError {
    /// A header record is not on the first line.
    UnexpectedHeader,
    /// The record is invalid.
    InvalidRecord(record::ParseError),
    /// A header record is invalid.
    InvalidHeader(header::TryFromRecordError),
    /// A reference sequence record is invalid.
    InvalidReferenceSequence(reference_sequence::TryFromRecordError),
    /// A reference read group record is invalid.
    InvalidReadGroup(read_group::TryFromRecordError),
    /// A program record is invalid.
    InvalidProgram(program::TryFromRecordError),
    /// A comment record is invalid.
    InvalidComment,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedHeader => f.write_str("unexpected @HD"),
            Self::InvalidRecord(e) => write!(f, "invalid record: {}", e),
            Self::InvalidHeader(e) => write!(f, "invalid header: {}", e),
            Self::InvalidReferenceSequence(e) => write!(f, "invalid reference sequence: {}", e),
            Self::InvalidReadGroup(e) => write!(f, "invalid read group: {}", e),
            Self::InvalidProgram(e) => write!(f, "invalid program: {}", e),
            Self::InvalidComment => f.write_str("invalid comment record"),
        }
    }
}

impl FromStr for Header {
    type Err = ParseError;

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
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut builder = Self::builder();

        for (i, line) in s.lines().enumerate() {
            let record: Record = line.parse().map_err(ParseError::InvalidRecord)?;

            builder = match record.kind() {
                record::Kind::Header => {
                    if i == 0 {
                        builder.set_header(
                            header::Header::try_from(record).map_err(ParseError::InvalidHeader)?,
                        )
                    } else {
                        return Err(ParseError::UnexpectedHeader);
                    }
                }
                record::Kind::ReferenceSequence => {
                    let reference_sequence = ReferenceSequence::try_from(record)
                        .map_err(ParseError::InvalidReferenceSequence)?;
                    builder.add_reference_sequence(reference_sequence)
                }
                record::Kind::ReadGroup => {
                    let read_group =
                        ReadGroup::try_from(record).map_err(ParseError::InvalidReadGroup)?;
                    builder.add_read_group(read_group)
                }
                record::Kind::Program => {
                    let program = Program::try_from(record).map_err(ParseError::InvalidProgram)?;
                    builder.add_program(program)
                }
                record::Kind::Comment => match record.value() {
                    record::Value::String(comment) => builder.add_comment(comment),
                    _ => return Err(ParseError::InvalidComment),
                },
            };
        }

        Ok(builder.build())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let header = Header::builder()
            .set_header(header::Header::new(String::from("1.6")))
            .add_reference_sequence(ReferenceSequence::new(String::from("sq0"), 8))
            .add_reference_sequence(ReferenceSequence::new(String::from("sq1"), 13))
            .add_read_group(ReadGroup::new(String::from("rg0")))
            .add_read_group(ReadGroup::new(String::from("rg1")))
            .add_program(Program::new(String::from("pg0")))
            .add_program(Program::new(String::from("pg1")))
            .add_comment("noodles")
            .add_comment("sam")
            .build();

        let actual = header.to_string();
        let expected = "\
@HD\tVN:1.6
@SQ\tSN:sq0\tLN:8
@SQ\tSN:sq1\tLN:13
@RG\tID:rg0
@RG\tID:rg1
@PG\tID:pg0
@PG\tID:pg1
@CO\tnoodles
@CO\tsam
";

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let s = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:sq0\tLN:1
@SQ\tSN:sq1\tLN:2
@RG\tID:rg0
@PG\tID:pg0\tPN:noodles
@CO\tnoodles_sam::header::tests::test_from_str
";

        let header: Header = s.parse()?;

        assert_eq!(header.reference_sequences().len(), 2);
        assert_eq!(header.read_groups().len(), 1);
        assert_eq!(header.programs().len(), 1);
        assert_eq!(header.comments.len(), 1);
        assert_eq!(
            &header.comments[0],
            "noodles_sam::header::tests::test_from_str"
        );

        Ok(())
    }

    #[test]
    fn test_from_str_with_empty_input() -> Result<(), ParseError> {
        let header: Header = "".parse()?;

        assert!(header.header().is_none());
        assert!(header.reference_sequences().is_empty());
        assert!(header.read_groups().is_empty());
        assert!(header.programs().is_empty());
        assert!(header.comments().is_empty());

        Ok(())
    }

    #[test]
    fn test_from_str_without_hd() -> Result<(), ParseError> {
        let header: Header = "@SQ\tSN:sq0\tLN:8\n".parse()?;
        assert!(header.header().is_none());
        assert_eq!(header.reference_sequences().len(), 1);
        Ok(())
    }

    #[test]
    fn test_from_str_with_multiple_hd() {
        let s = "\
@HD\tVN:1.6\tSO:coordinate
@HD\tVN:1.6\tSO:coordinate
";

        assert!(s.parse::<Header>().is_err());
    }
}
