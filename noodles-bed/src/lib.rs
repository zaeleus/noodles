//! **noodles-bed** handles the reading and writing of the BED (Browser Extensible Data) format.

pub mod feature;
pub mod fs;
pub mod io;
mod record;

pub use self::record::Record;
