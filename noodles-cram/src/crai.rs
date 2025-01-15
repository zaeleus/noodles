//! CRAM index.

#[cfg(feature = "async")]
pub mod r#async;

pub mod fs;
pub mod io;
pub mod record;

pub use self::record::Record;

#[deprecated(since = "0.74.0", note = "Use `crai::fs::read` instead.")]
pub use self::fs::read;

#[deprecated(since = "0.74.0", note = "Use `crai::fs::write` instead.")]
pub use self::fs::write;

#[deprecated(since = "0.74.0", note = "Use `crai::io::Reader` instead.")]
pub use self::io::Reader;

#[deprecated(since = "0.74.0", note = "Use `crai::io::Writer` instead.")]
pub use self::io::Writer;

#[cfg(feature = "async")]
#[deprecated(since = "0.74.0", note = "Use `crai::r#async::io::Reader` instead.")]
pub use self::r#async::io::Reader as AsyncReader;

#[cfg(feature = "async")]
#[deprecated(since = "0.74.0", note = "Use `crai::r#async::io::Writer` instead.")]
pub use self::r#async::io::Writer as AsyncWriter;

/// A CRAM index.
pub type Index = Vec<Record>;
