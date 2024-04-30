use std::io;

/// A buffered BGZF reader.
pub trait BufRead: io::BufRead + super::Read {}
