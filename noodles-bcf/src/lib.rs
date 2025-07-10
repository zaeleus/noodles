//! **noodles-bcf** handles the reading and writing of the BCF format.

#[cfg(feature = "async")]
pub mod r#async;

pub mod fs;
pub mod io;
pub mod record;

pub use self::record::Record;

#[cfg(feature = "async")]
#[deprecated(since = "0.76.0", note = "Use `bcf::r#async::io::Reader` instead.")]
pub use self::r#async::io::Reader as AsyncReader;

#[cfg(feature = "async")]
#[deprecated(since = "0.76.0", note = "Use `bcf::r#async::io::Writer` instead.")]
pub use self::r#async::io::Writer as AsyncWriter;
