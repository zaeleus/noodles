//! SAM header record.

pub mod kind;
pub mod value;

use bstr::BString;

pub use self::kind::Kind;
use self::value::{
    map::{self, Program, ReadGroup, ReferenceSequence},
    Map,
};

/// A SAM header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Record {
    /// A header (`HD`) record.
    Header(Map<map::Header>),
    /// A reference sequence (`SQ`) record.
    ReferenceSequence(Vec<u8>, Map<ReferenceSequence>),
    /// A read group (`RG`) record.
    ReadGroup(BString, Map<ReadGroup>),
    /// A program (`PG`) record.
    Program(BString, Map<Program>),
    /// A comment (`CO`) record.
    Comment(Vec<u8>),
}
