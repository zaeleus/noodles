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
//! # #[cfg(feature = "records")]
//! # fn main() -> std::io::Result<()> {
//! # use std::fs::File;
//! use noodles_bam as bam;
//!
//! let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
//! let header = reader.read_header()?;
//!
//! for result in reader.records() {
//!     let record = result?;
//!     // ...
//! }
//! # Ok(())
//! # }
//! # #[cfg(not(feature = "records"))]
//! # fn main() {}
//! ```
//!
//! ## Query records
//!
//! Querying allows filtering records by region. It requires an associated BAM index (BAI).
//!
//! ```no_run
//! # #[cfg(feature = "records")]
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! use noodles_bam as bam;
//!
//! let mut reader = bam::io::indexed_reader::Builder::default().build_from_path("sample.bam")?;
//! let header = reader.read_header()?;
//!
//! let region = "sq0:5-8".parse()?;
//! let query = reader.query(&header, &region)?;
//!
//! for result in query.records() {
//!     let record = result?;
//!     // ...
//! }
//! # Ok(())
//! # }
//! # #[cfg(not(feature = "records"))]
//! # fn main() {}
//! ```

#[cfg(feature = "async")]
pub mod r#async;

pub mod bai;
#[cfg(feature = "records")]
pub mod fs;
#[cfg(feature = "records")]
pub mod io;
#[cfg(feature = "records")]
pub mod record;
#[cfg(feature = "records")]
mod record_ref;

#[cfg(feature = "records")]
pub use self::{record::Record, record_ref::RecordRef};
