//! Async GFF.

pub mod io;

#[deprecated(
    since = "0.33.0",
    note = "Use `noodles_gff::r#async::io::Reader` instead."
)]
pub use self::io::Reader;
