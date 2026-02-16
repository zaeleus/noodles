//! Indexed CRAM reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, Read, Seek};

use noodles_core::Region;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::{
    Reader,
    reader::{Container, Query, Records},
};
use crate::{FileDefinition, crai};

/// An indexed CRAM reader.
pub struct IndexedReader<R> {
    inner: Reader<R>,
    index: crai::Index,
}

impl<R> IndexedReader<R>
where
    R: Read,
{
    /// Creates an indexed CRAM reader.
    pub fn new(inner: R, index: crai::Index) -> Self {
        Self {
            inner: Reader::new(inner),
            index,
        }
    }

    /// Returns a reference to the underlying reader.
    pub fn get_ref(&self) -> &R {
        self.inner.get_ref()
    }

    /// Returns a mutable reference to the underlying reader.
    pub fn get_mut(&mut self) -> &mut R {
        self.inner.get_mut()
    }

    /// Unwraps and returns the underlying reader.
    pub fn into_inner(self) -> R {
        self.inner.into_inner()
    }

    /// Returns the reference sequence repository.
    pub fn reference_sequence_repository(&self) -> &fasta::Repository {
        self.inner.reference_sequence_repository()
    }

    /// Reads the CRAM file definition.
    pub fn read_file_definition(&mut self) -> io::Result<FileDefinition> {
        self.inner.read_file_definition()
    }

    /// Reads the SAM header.
    pub fn read_file_header(&mut self) -> io::Result<sam::Header> {
        self.inner.read_file_header()
    }

    /// Reads the SAM header.
    pub fn read_header(&mut self) -> io::Result<sam::Header> {
        self.inner.read_header()
    }

    /// Reads a container.
    pub fn read_container(&mut self, container: &mut Container) -> io::Result<usize> {
        self.inner.read_container(container)
    }

    /// Returns a iterator over records starting from the current stream position.
    pub fn records<'r, 'h: 'r>(&'r mut self, header: &'h sam::Header) -> Records<'r, 'h, R> {
        self.inner.records(header)
    }

    /// Returns the associated index.
    pub fn index(&self) -> &crai::Index {
        &self.index
    }

    /// Calls `f` for each record, avoiding intermediate `RecordBuf` conversion.
    ///
    /// This delegates to [`Reader::for_each_record`] and is useful when writing CRAM records
    /// directly to another format without first converting to `RecordBuf`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_cram as cram;
    ///
    /// let mut reader = cram::io::indexed_reader::Builder::default()
    ///     .build_from_path("sample.cram")?;
    /// let header = reader.read_header()?;
    ///
    /// reader.for_each_record(&header, |record| {
    ///     // record is &dyn sam::alignment::Record
    ///     Ok(())
    /// })?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn for_each_record<F>(&mut self, header: &sam::Header, f: F) -> io::Result<()>
    where
        F: FnMut(&dyn sam::alignment::Record) -> io::Result<()>,
    {
        self.inner.for_each_record(header, f)
    }
}

impl<R> IndexedReader<R>
where
    R: Read + Seek,
{
    /// Returns an iterator over records that intersect the given region.
    ///
    /// The index owned by this reader is used to find the relevant containers. This is equivalent
    /// to calling [`Reader::query`] with the associated index.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_cram as cram;
    ///
    /// let mut reader = cram::io::indexed_reader::Builder::default()
    ///     .build_from_path("sample.cram")?;
    /// let header = reader.read_header()?;
    ///
    /// let region = "sq0:8-13".parse()?;
    ///
    /// for result in reader.query(&header, &region)? {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        region: &Region,
    ) -> io::Result<Query<'r, 'h, 'r, R>> {
        self.inner.query(header, &self.index, region)
    }

    /// Returns an iterator over unmapped records.
    ///
    /// The index owned by this reader is used to seek to the unmapped records. This is equivalent
    /// to calling [`Reader::query_unmapped`] with the associated index.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_cram as cram;
    ///
    /// let mut reader = cram::io::indexed_reader::Builder::default()
    ///     .build_from_path("sample.cram")?;
    /// let header = reader.read_header()?;
    ///
    /// for result in reader.query_unmapped(&header)? {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn query_unmapped<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
    ) -> io::Result<impl Iterator<Item = io::Result<sam::alignment::RecordBuf>> + use<'r, 'h, R>>
    {
        self.inner.query_unmapped(header, &self.index)
    }
}
