//! Indexed alignment reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, Read, Seek};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_cram as cram;
use noodles_sam::{self as sam, alignment::Record};

/// An indexed alignment reader.
pub enum IndexedReader<R> {
    /// SAM.
    Sam(sam::io::IndexedReader<R>),
    /// BAM.
    Bam(bam::io::IndexedReader<bgzf::io::Reader<R>>),
    /// CRAM.
    Cram(cram::io::IndexedReader<R>),
}

impl<R> IndexedReader<R>
where
    R: Read,
{
    /// Reads the SAM header.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_util::alignment;
    ///
    /// let mut reader = alignment::io::indexed_reader::Builder::default()
    ///     .build_from_path("sample.bam")?;
    ///
    /// let _header = reader.read_header()?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<sam::Header> {
        match self {
            Self::Sam(reader) => reader.read_header(),
            Self::Bam(reader) => reader.read_header(),
            Self::Cram(reader) => reader.read_header(),
        }
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_util::alignment;
    ///
    /// let mut reader = alignment::io::indexed_reader::Builder::default()
    ///     .build_from_path("sample.bam")?;
    ///
    /// let header = reader.read_header()?;
    ///
    /// for result in reader.records(&header) {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
    ) -> impl Iterator<Item = io::Result<Box<dyn Record>>> + use<'r, 'h, R> {
        let records: Box<dyn Iterator<Item = io::Result<Box<dyn Record>>>> = match self {
            Self::Sam(reader) => Box::new(
                reader
                    .records()
                    .map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
            ),
            Self::Bam(reader) => Box::new(
                reader
                    .records()
                    .map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
            ),
            Self::Cram(reader) => Box::new(
                reader
                    .records(header)
                    .map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
            ),
        };

        records
    }
}

impl<R> IndexedReader<R>
where
    R: Read + Seek,
{
    /// Returns an iterator over records that intersects the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_util::alignment;
    ///
    /// let mut reader = alignment::io::indexed_reader::Builder::default()
    ///     .build_from_path("sample.bam")?;
    ///
    /// let header = reader.read_header()?;
    ///
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        region: &Region,
    ) -> io::Result<impl Iterator<Item = io::Result<Box<dyn Record>>> + use<'r, 'h, R>> {
        let records: Box<dyn Iterator<Item = io::Result<Box<dyn Record>>>> = match self {
            Self::Sam(reader) => {
                let query = reader.query(header, region)?;

                Box::new(
                    query.map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
                )
            }
            Self::Bam(reader) => {
                let query = reader.query(header, region)?;

                Box::new(
                    query
                        .records()
                        .map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
                )
            }
            Self::Cram(reader) => {
                let query = reader.query(header, region)?;

                Box::new(
                    query.map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
                )
            }
        };

        Ok(records)
    }
}
