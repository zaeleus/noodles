#![warn(missing_docs)]

//! **noodles-vcf** handles the reading and writing of the VCF format.
//!
//! # Examples
//!
//! ## Read all records from a file
//!
//! ```no_run
//! # use std::{fs::File, io::{self, BufReader}};
//! use noodles_vcf as vcf;
//!
//! let mut reader = File::open("sample.vcf").map(BufReader::new).map(vcf::Reader::new)?;
//! reader.read_header()?;
//!
//! for result in reader.records() {
//!     let record = result?;
//!     println!("{:?}", record);
//! }
//! # Ok::<(), io::Error>(())
//! ```

#[cfg(feature = "async")]
mod r#async;

pub mod header;
mod reader;
pub mod record;
mod writer;

pub use self::{header::Header, reader::Reader, record::Record, writer::Writer};

#[cfg(feature = "async")]
pub use self::r#async::Reader as AsyncReader;
