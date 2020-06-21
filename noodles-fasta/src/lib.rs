//! **noodles-fasta** handles and reading and writing of the FASTA format.
//!
//! FASTA is a text format with no formal specification and only has de facto rules. It typically
//! consists of a list of records, each with a definition on the first line and a sequence in the
//! following lines.
//!
//! The definition starts with a `>` (greater than) character, and directly after it is the
//! reference sequence name. Optionally, whitespace may be used a delimiter for an extra
//! description or metadata of the sequence. For example,
//!
//! ```text
//!  reference sequence name
//!  | |
//! >sq0 LN:13
//!      |   |
//!      description
//! ```
//!
//! The sequence is effectively a byte array of characters representing a base. It is typically
//! hard wrapped at an arbitrary width. For example, the following makes up the sequence
//! `ACGTNACTGG`.
//!
//! ```text
//! ACGT
//! NACT
//! GG
//! ```

pub mod fai;
mod reader;
pub mod record;

pub use self::{reader::Reader, record::Record};
