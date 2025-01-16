//! Async gzip index.

pub mod fs;
pub mod io;

#[deprecated(since = "0.35.0", note = "Use `gzi::r#async::fs::read` instead.")]
pub use self::fs::read;

#[deprecated(since = "0.35.0", note = "Use `gzi::r#async::io::Reader` instead.")]
pub use self::io::Reader;
