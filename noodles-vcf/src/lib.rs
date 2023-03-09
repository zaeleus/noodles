#![warn(missing_docs)]

//! **noodles-vcf** handles the reading and writing of the VCF format.
//!
//! # Examples
//!
//! ## Read all records from a file
//!
//! ```no_run
//! # use std::{fs::File, io::BufReader};
//! use noodles_vcf as vcf;
//!
//! let mut reader = File::open("sample.vcf").map(BufReader::new).map(vcf::Reader::new)?;
//! let header = reader.read_header()?.parse()?;
//!
//! for result in reader.records(&header) {
//!     let record = result?;
//!     println!("{:?}", record);
//! }
//! # Ok::<_, Box<dyn std::error::Error>>(())
//! ```

#[cfg(feature = "async")]
mod r#async;

pub mod header;
pub mod indexed_reader;
pub mod reader;
pub mod record;
mod variant_reader;
mod writer;

pub use self::{
    header::Header, indexed_reader::IndexedReader, reader::Reader, record::Record,
    variant_reader::VariantReader, writer::Writer,
};

#[cfg(feature = "async")]
pub use self::r#async::{Reader as AsyncReader, Writer as AsyncWriter};
