//! **noodles-bcf** handles the reading and writing of the BCF format.

pub mod header;
mod reader;

pub use self::reader::Reader;
