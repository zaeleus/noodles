//! Indexed CRAM reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, Read, Seek};

use noodles_core::Region;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::{
    crai,
    reader::{Query, Records},
    DataContainer, FileDefinition, Reader,
};

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

    /// Reads the CRAM file definition.
    pub fn read_file_definition(&mut self) -> io::Result<FileDefinition> {
        self.inner.read_file_definition()
    }

    /// Reads the raw SAM header.
    pub fn read_file_header(&mut self) -> io::Result<sam::Header> {
        self.inner.read_file_header()
    }

    /// Reads a data container.
    pub fn read_data_container(&mut self) -> io::Result<Option<DataContainer>> {
        self.inner.read_data_container()
    }

    /// Returns a iterator over records starting from the current stream position.
    pub fn records<'a>(
        &'a mut self,
        reference_sequence_repository: &'a fasta::Repository,
        header: &'a sam::Header,
    ) -> Records<'a, R> {
        self.inner.records(reference_sequence_repository, header)
    }

    /// Returns the associated index.
    pub fn index(&self) -> &crai::Index {
        &self.index
    }
}

impl<R> IndexedReader<R>
where
    R: Read + Seek,
{
    /// Returns an iterator over records that intersects the given region.
    pub fn query<'a>(
        &'a mut self,
        reference_sequence_repository: &'a fasta::Repository,
        header: &'a sam::Header,
        region: &Region,
    ) -> io::Result<Query<'_, R>> {
        self.inner
            .query(reference_sequence_repository, header, &self.index, region)
    }
}
