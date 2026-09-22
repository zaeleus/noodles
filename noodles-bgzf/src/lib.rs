//! **noodles-bgzf** handles the reading and writing of the blocked gzip format (BGZF).
//!
//! While the gzip format is typically a single stream, a BGZF is the concatenation of many gzip
//! streams. Each stream is called a block, with its uncompressed data size being constrained to
//! less than 64 KiB. This multistream gzip allows random access using [`virtual positions`].
//!
//! noodles-bgzf abstracts away the concept of blocks, implementing [`std::io::Read`] for the
//! reader and [`std::io::Write`] for the writer.
//!
//! [`virtual positions`]: VirtualPosition
//!
//! # Examples
//!
//! ## Read an entire BGZF file
//!
//! ```no_run
//! # use std::{fs::File, io::{self, Read}};
//! use noodles_bgzf as bgzf;
//! let mut reader = File::open("data.gz").map(bgzf::io::Reader::new)?;
//! let mut data = Vec::new();
//! reader.read_to_end(&mut data)?;
//! # Ok::<(), io::Error>(())
//! ```
//!
//! ## Write a BGZF file
//!
//! ```no_run
//! # use std::{fs::File, io::{self, Write}};
//! use noodles_bgzf as bgzf;
//! let mut writer = File::create("data.gz").map(bgzf::io::Writer::new)?;
//! writer.write_all(b"noodles-bgzf")?;
//! # Ok::<(), io::Error>(())
//! ```

#[cfg(feature = "async")]
pub mod r#async;

pub(crate) mod deflate;
pub mod fs;
mod gz;
pub mod gzi;
pub mod io;
pub mod virtual_position;

pub use self::virtual_position::VirtualPosition;

// XLEN (2)
const GZIP_XLEN_SIZE: usize = 2;

// SI1 (1) + SI2 (1) + SLEN (2) + BSIZE (2)
const BGZF_XLEN: usize = 6;

// § 4.1 The BGZF compression format (2021-06-03): "Thus while `ISIZE` is stored as a `uint32_t` as
// per the gzip format, in BGZF it is limited to the range [0, 65536]."
const BGZF_MAX_ISIZE: usize = 1 << 16;

pub(crate) const BGZF_HEADER_SIZE: usize = gz::HEADER_SIZE + GZIP_XLEN_SIZE + BGZF_XLEN;
