//! SAM header and records.
//!
//! A SAM header is a list of header records. There are 5 record types: [header] (`HD`), [reference
//! sequence] (`SQ`), [read group] (`RG`), [program] (`PG`), and comment (`CO`).
//!
//! Each record is effectively a map. It defines key-value pairs associated with that record type.
//!
//! All records are optional, which means an empty header is considered a valid SAM header.
//!
//! If there is a header record, it must appear as the first record.
//!
//! Reference sequence, read group, program, and comment records are lists of records of the same
//! type. Reference sequences must be ordered; whereas read groups, programs, and comments can be
//! unordered. (`sam::Header` defines them to be ordered.)
//!
//! [header]: `record::value::map::Header`
//! [reference sequence]: `ReferenceSequence`
//! [read group]: `ReadGroup`
//! [program]: `Program`
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
//! use std::num::NonZeroUsize;
//!
//! use noodles_sam::{
//!     self as sam,
//!     header::record::value::{map::ReferenceSequence, Map},
//! };
//!
//! let header = sam::Header::builder()
//!     .set_header(Default::default())
//!     .add_reference_sequence(
//!         "sq0",
//!         Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8)?),
//!     )
//!     .add_reference_sequence(
//!         "sq1",
//!         Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?),
//!     )
//!     .build();
//!
//! assert!(header.header().is_some());
//! assert_eq!(header.reference_sequences().len(), 2);
//! assert!(header.read_groups().is_empty());
//! assert!(header.programs().is_empty());
//! assert!(header.comments().is_empty());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod builder;
mod parser;
pub mod record;

pub use self::{
    builder::Builder,
    parser::{ParseError, Parser},
    record::Record,
};

use std::str::{self, FromStr};

use bstr::BString;
use indexmap::IndexMap;

use self::record::value::{
    map::{self, Program, ReadGroup, ReferenceSequence},
    Map,
};

/// A reference sequence dictionary.
pub type ReferenceSequences = IndexMap<BString, Map<ReferenceSequence>>;

/// An ordered map of read groups.
pub type ReadGroups = IndexMap<BString, Map<ReadGroup>>;

/// An ordered map of programs.
pub type Programs = IndexMap<BString, Map<Program>>;

/// A SAM header.
///
/// Records are grouped by their types: header, reference sequence, read group, program, and
/// comment.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Header {
    header: Option<Map<map::Header>>,
    reference_sequences: ReferenceSequences,
    read_groups: ReadGroups,
    programs: Programs,
    comments: Vec<Vec<u8>>,
}

impl Header {
    /// Returns a builder to create a SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let builder = sam::Header::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns the SAM header header if it is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     self as sam,
    ///     header::record::value::{map::Header, Map},
    /// };
    ///
    /// let header = sam::Header::default();
    /// assert!(header.header().is_none());
    ///
    /// let header = sam::Header::builder()
    ///     .set_header(Map::<Header>::default())
    ///     .build();
    ///
    /// assert!(header.header().is_some());
    /// ```
    pub fn header(&self) -> Option<&Map<map::Header>> {
        self.header.as_ref()
    }

    /// Returns a mutable reference to the SAM header header if it is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::record::value::{map, Map}};
    ///
    /// let mut header = sam::Header::default();
    /// assert!(header.header().is_none());
    ///
    /// let hd = Map::<map::Header>::default();
    /// *header.header_mut() = Some(hd.clone());
    /// assert_eq!(header.header(), Some(&hd));
    /// ```
    pub fn header_mut(&mut self) -> &mut Option<Map<map::Header>> {
        &mut self.header
    }

    /// Returns the SAM header reference sequences.
    ///
    /// This is also called the reference sequence dictionary.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::{
    ///     self as sam,
    ///     header::record::value::{map::ReferenceSequence, Map},
    /// };
    ///
    /// let header = sam::Header::builder()
    ///     .add_reference_sequence(
    ///         "sq0",
    ///         Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?)
    ///     )
    ///     .build();
    ///
    /// let reference_sequences = header.reference_sequences();
    /// assert_eq!(reference_sequences.len(), 1);
    /// assert!(reference_sequences.contains_key(&b"sq0"[..]));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
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
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::{
    ///     self as sam,
    ///     header::record::value::{map::ReferenceSequence, Map},
    /// };
    ///
    /// let mut header = sam::Header::default();
    ///
    /// header.reference_sequences_mut().insert(
    ///     String::from("sq0").into(),
    ///     Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?),
    /// );
    ///
    /// let reference_sequences = header.reference_sequences();
    /// assert_eq!(reference_sequences.len(), 1);
    /// assert!(reference_sequences.contains_key(&b"sq0"[..]));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_sequences_mut(&mut self) -> &mut ReferenceSequences {
        &mut self.reference_sequences
    }

    /// Returns the SAM header read groups.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     self as sam,
    ///     header::record::value::{map::ReadGroup, Map},
    /// };
    ///
    /// let header = sam::Header::builder()
    ///     .add_read_group("rg0", Map::<ReadGroup>::default())
    ///     .build();
    ///
    /// let read_groups = header.read_groups();
    /// assert_eq!(read_groups.len(), 1);
    /// assert!(read_groups.contains_key(&b"rg0"[..]));
    /// ```
    pub fn read_groups(&self) -> &ReadGroups {
        &self.read_groups
    }

    /// Returns a mutable reference to the SAM header read groups.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     self as sam,
    ///     header::record::value::{map::ReadGroup, Map},
    /// };
    ///
    /// let mut header = sam::Header::default();
    /// assert!(header.read_groups().is_empty());
    ///
    /// let read_group = Map::<ReadGroup>::default();
    /// header.read_groups_mut().insert(String::from("rg0").into(), read_group);
    ///
    /// let read_groups = header.read_groups();
    /// assert_eq!(read_groups.len(), 1);
    /// assert!(read_groups.contains_key(&b"rg0"[..]));
    /// ```
    pub fn read_groups_mut(&mut self) -> &mut ReadGroups {
        &mut self.read_groups
    }

    /// Returns the SAM header programs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::record::value::{map::Program, Map}};
    ///
    /// let program = Map::<Program>::default();
    /// let header = sam::Header::builder().add_program("noodles-sam", program).build();
    ///
    /// let programs = header.programs();
    /// assert_eq!(programs.len(), 1);
    /// assert!(programs.contains_key(&b"noodles-sam"[..]));
    /// ```
    pub fn programs(&self) -> &Programs {
        &self.programs
    }

    /// Returns a mutable reference to the SAM header programs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::record::value::{map::Program, Map}};
    ///
    /// let mut header = sam::Header::default();
    ///
    /// let program = Map::<Program>::default();
    /// header.programs_mut().insert(String::from("noodles-sam").into(), program);
    ///
    /// let programs = header.programs();
    /// assert_eq!(programs.len(), 1);
    /// assert!(programs.contains_key(&b"noodles-sam"[..]));
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
    /// assert_eq!(header.comments(), [Vec::from("noodles-sam")]);
    /// ```
    pub fn comments(&self) -> &[Vec<u8>] {
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
    /// let mut header = sam::Header::default();
    /// header.comments_mut().push(Vec::from("noodles-sam"));
    /// assert_eq!(header.comments(), [Vec::from("noodles-sam")]);
    /// ```
    pub fn comments_mut(&mut self) -> &mut Vec<Vec<u8>> {
        &mut self.comments
    }

    /// Adds a comment.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let mut header = sam::Header::default();
    /// header.add_comment("noodles-sam");
    /// assert_eq!(header.comments(), [Vec::from("noodles-sam")]);
    /// ```
    pub fn add_comment<C>(&mut self, comment: C)
    where
        C: Into<Vec<u8>>,
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

    /// Removes all records from the header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let mut header = sam::Header::builder().add_comment("ndls").build();
    /// assert!(!header.is_empty());
    ///
    /// header.clear();
    /// assert!(header.is_empty());
    /// ```
    pub fn clear(&mut self) {
        self.header.take();
        self.reference_sequences.clear();
        self.read_groups.clear();
        self.programs.clear();
        self.comments.clear();
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
        parser::parse(s)
    }
}
