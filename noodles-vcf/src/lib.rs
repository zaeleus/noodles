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
//! let mut reader = vcf::reader::Builder::default().build_from_path("sample.vcf")?;
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
pub mod indexed_reader;
mod indexer;
pub mod io;
pub mod lazy;
pub mod reader;
pub mod record;
mod variant_reader;
mod variant_writer;
pub mod writer;

pub use self::{
    header::Header, indexed_reader::IndexedReader, indexer::index, reader::Reader, record::Record,
    variant_reader::VariantReader, variant_writer::VariantWriter, writer::Writer,
};

#[cfg(feature = "async")]
pub use self::r#async::{Reader as AsyncReader, Writer as AsyncWriter};
