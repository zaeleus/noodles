//! FASTA index (FAI) and fields.

#[cfg(feature = "async")]
pub mod r#async;

pub mod fs;
mod index;
pub mod io;
mod record;

pub use self::{index::Index, record::Record};

#[deprecated(since = "0.47.0", note = "Use `fai::fs::read` instead.")]
pub use self::fs::read;

#[deprecated(since = "0.44.0", note = "Use `fai::io::Reader` instead.")]
pub use self::io::Reader;

#[deprecated(since = "0.44.0", note = "Use `fai::io::Writer` instead.")]
pub use self::io::Writer;

#[cfg(feature = "async")]
#[deprecated(since = "0.44.0", note = "Use `fai::r#async::io::Reader` instead.")]
pub use self::r#async::io::Reader as AsyncReader;
