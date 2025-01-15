//! Async BAI.

pub mod fs;
pub mod io;

#[deprecated(since = "0.73.0", note = "Use `bai::r#async::fs::read` instead.")]
pub use self::fs::read;

#[deprecated(since = "0.73.0", note = "Use `bai::r#async::fs::write` instead.")]
pub use self::fs::write;

#[deprecated(since = "0.68.0", note = "Use `bai::r#async::io::Reader` instead.")]
pub use self::io::Reader;

#[deprecated(since = "0.68.0", note = "Use `bai::r#async::io::Writer` instead.")]
pub use self::io::Writer;
