//! BCF header.

pub mod string_maps;

pub use self::string_maps::StringMaps;

#[deprecated(
    since = "0.11.0",
    note = "Use `noodles_bcf::header::string_maps::StringMap` instead."
)]
pub use self::string_maps::StringMap;
