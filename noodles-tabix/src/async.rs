//! Async tabix.

pub mod fs;
pub mod io;

#[deprecated(since = "0.48.0", note = "Use `tabix::r#async::fs::read` instead.")]
pub use self::fs::read;

#[deprecated(since = "0.48.0", note = "Use `tabix::r#async::fs::write` instead.")]
pub use self::fs::write;

#[deprecated(
    since = "0.45.0",
    note = "Use `noodles_tabix::r#async::io::Reader` instead."
)]
pub use self::io::Reader;

#[deprecated(
    since = "0.45.0",
    note = "Use `noodles_tabix::r#async::io::Reader` instead."
)]
pub use self::io::Writer;
