//! Indexed BAM reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_sam::{alignment::Record, header::ReferenceSequences};

use crate::reader::UnmappedRecords;

use super::{
    bai, lazy,
    reader::{LazyRecords, Query, Records},
    Reader,
};

/// An indexed BAM reader.
pub struct IndexedReader<R> {
    inner: Reader<R>,
    index: bai::Index,
}

impl<R> IndexedReader<R>
where
    R: Read,
{
    /// Returns a reference to the underlying reader.
    pub fn get_ref(&self) -> &R {
        self.inner.get_ref()
    }

    /// Returns a mutable reference to the underlying reader.
    pub fn get_mut(&mut self) -> &mut R {
        self.inner.get_mut()
    }

    /// Returns the underlying reader.
    pub fn into_inner(self) -> R {
        self.inner.into_inner()
    }

    /// Reads the raw SAM header.
    pub fn read_header(&mut self) -> io::Result<String> {
        self.inner.read_header()
    }

    /// Reads the binary reference sequences after the SAM header.
    pub fn read_reference_sequences(&mut self) -> io::Result<ReferenceSequences> {
        self.inner.read_reference_sequences()
    }

    /// Reads a single record.
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        self.inner.read_record(record)
    }

    /// Reads a single record without eagerly decoding its fields.
    pub fn read_lazy_record(&mut self, record: &mut lazy::Record) -> io::Result<usize> {
        self.inner.read_lazy_record(record)
    }

    /// Returns an iterator over records starting from the current stream position.
    pub fn records(&mut self) -> Records<'_, R> {
        self.inner.records()
    }

    /// Returns an iterator over lazy records.
    pub fn lazy_records(&mut self) -> LazyRecords<'_, R> {
        self.inner.lazy_records()
    }
}

impl<R> IndexedReader<bgzf::Reader<R>>
where
    R: Read,
{
    /// Creates an indexed BAM reader.
    pub fn new(inner: R, index: bai::Index) -> Self {
        Self {
            inner: Reader::new(inner),
            index,
        }
    }
}

impl<R> IndexedReader<bgzf::Reader<R>>
where
    R: Read + Seek,
{
    /// Returns an iterator over records that intersect the given region.
    pub fn query(
        &mut self,
        reference_sequences: &ReferenceSequences,
        region: &Region,
    ) -> io::Result<Query<'_, R>> {
        self.inner.query(reference_sequences, &self.index, region)
    }

    /// Returns an iterator of unmapped records after querying for the unmapped region.
    pub fn query_unmapped(&mut self) -> io::Result<UnmappedRecords<'_, R>> {
        self.inner.query_unmapped(&self.index)
    }
}
