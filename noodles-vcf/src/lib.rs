#![warn(missing_docs)]

//! **noodles-vcf** handles the reading and writing of the VCF format.
//!
//! # Examples
//!
//! ## Read all records from a file
//!
//! ```no_run
//! use noodles_vcf as vcf;
//!
//! let mut reader = vcf::io::reader::Builder::default().build_from_path("sample.vcf")?;
//! let header = reader.read_header()?;
//!
//! for result in reader.records() {
//!     let record = result?;
//!     // ...
//! }
//! # Ok::<_, std::io::Error>(())
//! ```

#[cfg(feature = "async")]
pub mod r#async;

pub mod fs;
pub mod header;
pub mod io;
pub mod record;
pub mod variant;

pub use self::{header::Header, record::Record};

#[deprecated(since = "0.73.0", note = "Use `vcf::fs::index` instead.")]
pub use self::fs::index;

#[cfg(feature = "async")]
pub use self::r#async::io::{Reader as AsyncReader, Writer as AsyncWriter};
