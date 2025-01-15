//! Async CSI.

pub mod fs;
pub mod io;

#[deprecated(since = "0.42.0", note = "Use `csi::r#async::fs::read` instead.")]
pub use self::fs::read;

#[deprecated(since = "0.42.0", note = "Use `csi::r#async::fs::write` instead.")]
pub use self::fs::write;

#[deprecated(
    since = "0.39.0",
    note = "Use `noodles_csi::r#async::io::Reader` instead."
)]
pub use self::io::Reader;

#[deprecated(
    since = "0.39.0",
    note = "Use `noodles_csi::r#async::io::Writer` instead."
)]
pub use self::io::Writer;
