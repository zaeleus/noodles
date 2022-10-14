//! FASTA I/O.

use std::io::{BufRead, Seek};

/// A reader that is both buffered and seekable.
pub trait BufReadSeek: BufRead + Seek {}

impl<T> BufReadSeek for T where T: BufRead + Seek {}
