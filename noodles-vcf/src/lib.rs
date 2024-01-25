#![allow(deprecated)]
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
//! for result in reader.records(&header) {
//!     let record = result?;
//!     println!("{:?}", record);
//! }
//! # Ok::<_, std::io::Error>(())
//! ```

#[cfg(feature = "async")]
mod r#async;

pub mod header;
mod indexer;
pub mod io;
pub mod lazy;
pub mod record;
pub mod variant;

pub use self::{header::Header, indexer::index, record::Record};

#[cfg(feature = "async")]
pub use self::r#async::{Reader as AsyncReader, Writer as AsyncWriter};
