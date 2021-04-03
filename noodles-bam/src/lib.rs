#![deny(missing_docs)]

//! **noodles-bam** handles the reading and writing of the BAM (Binary Alignment/Map) file format.
//!
//! The BAM format contains the same information as SAM (Sequence Alignment/Map), namely a SAM
//! header and a list of records.
//!
//! # Examples
//!
//! ## Read all records
//!
//! ```no_run
//! # use std::{fs::File, io};
//! use noodles_bam as bam;
//!
//! let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
//! reader.read_header()?;
//! reader.read_reference_sequences()?;
//!
//! for result in reader.records() {
//!     let record = result?;
//!     println!("{:?}", record);
//! }
//! # Ok::<(), io::Error>(())
//! ```
//!
//! ## Query records
//!
//! Querying allows filtering records by region. It requires an associated BAM index (BAI).
//!
//! ```no_run
//! # use std::fs::File;
//! use noodles_bam::{self as bam, bai};
//! use noodles_core::Region;
//! use noodles_sam as sam;
//!
//! let mut reader = File::open("sample.bam").map(bam::Reader::new)?;
//! let header: sam::Header = reader.read_header()?.parse()?;
//!
//! let reference_sequences = header.reference_sequences();
//! let index = bai::read("sample.bam.bai")?;
//! let region = Region::mapped("sq0", 17711, 28657);
//! let query = reader.query(&reference_sequences, &index, &region)?;
//!
//! for result in query {
//!     let record = result?;
//!     println!("{:?}", record);
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub mod bai;
pub mod reader;
pub mod record;
mod writer;

pub use self::{reader::Reader, record::Record, writer::Writer};

static MAGIC_NUMBER: &[u8] = b"BAM\x01";
