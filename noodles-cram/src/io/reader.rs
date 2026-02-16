//! CRAM reader and record iterator.

mod builder;
pub(crate) mod collections;
pub(crate) mod container;
pub mod header;
pub(crate) mod num;
mod query;
mod records;

use std::io::{self, Read, Seek, SeekFrom};

use noodles_core::Region;
use noodles_fasta as fasta;
use noodles_sam as sam;

use self::container::read_container;
pub use self::{builder::Builder, container::Container, query::Query, records::Records};
use crate::{FileDefinition, crai, file_definition::Version};

/// A CRAM reader.
///
/// The CRAM format is comprised of four main parts: 1) a file definition, 2) a file header, 3) a
/// list of containers, and 4) an end-of-file (EOF) container.
///
/// # Examples
///
/// ```no_run
/// # use std::{fs::File, io};
/// use noodles_cram as cram;
/// use noodles_fasta as fasta;
///
/// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
/// let header = reader.read_header()?;
///
/// for result in reader.records(&header) {
///     let record = result?;
///     // ...
/// }
///
/// # Ok::<_, io::Error>(())
/// ```
pub struct Reader<R> {
    inner: R,
    reference_sequence_repository: fasta::Repository,
    version: Version,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let reader = cram::io::Reader::new(io::empty());
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let mut reader = cram::io::Reader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Unwraps and returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let reader = cram::io::Reader::new(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a CRAM reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(inner: R) -> Self {
        Builder::default().build_from_reader(inner)
    }

    pub(crate) fn reference_sequence_repository(&self) -> &fasta::Repository {
        &self.reference_sequence_repository
    }

    /// Returns a CRAM header reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io::Read};
    /// use noodles_cram as cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number()?;
    /// header_reader.read_format_version()?;
    /// header_reader.read_file_id()?;
    ///
    /// let mut container_reader = header_reader.container_reader()?;
    ///
    /// let _raw_header = {
    ///     let mut raw_sam_header_reader = container_reader.raw_sam_header_reader()?;
    ///     let mut raw_header = String::new();
    ///     raw_sam_header_reader.read_to_string(&mut raw_header)?;
    ///     raw_sam_header_reader.discard_to_end()?;
    ///     raw_header
    /// };
    ///
    /// container_reader.discard_to_end()?;
    /// Ok::<_, std::io::Error>(())
    /// ```
    pub fn header_reader(&mut self) -> header::Reader<&mut R> {
        header::Reader::new(&mut self.inner)
    }

    /// Reads the CRAM file definition.
    ///
    /// The CRAM magic number is also checked.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// let file_definition = reader.read_file_definition()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_file_definition(&mut self) -> io::Result<FileDefinition> {
        let file_definition = header::read_file_definition(&mut self.inner)?;
        self.version = file_definition.version();
        self.version.validate()?;
        Ok(file_definition)
    }

    /// Reads the SAM header.
    ///
    /// The position of the stream is expected to be at the CRAM header container, i.e., directly
    /// after the file definition.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// reader.read_file_definition()?;
    ///
    /// let header = reader.read_file_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_file_header(&mut self) -> io::Result<sam::Header> {
        header::read_file_header(&mut self.inner, self.version)
    }

    /// Reads the SAM header.
    ///
    /// This verifies the CRAM magic number, discards the file definition, and reads and parses the
    /// file header as a SAM header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// let header = reader.read_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<sam::Header> {
        self.read_file_definition()?;
        self.read_file_header()
    }

    /// Reads a container.
    ///
    /// This returns `None` if the container header is the EOF container header, which signals the
    /// end of the stream.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram::{self as cram, io::reader::Container};
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// reader.read_header()?;
    ///
    /// let mut container = Container::default();
    ///
    /// while reader.read_container(&mut container)? != 0 {
    ///     // ...
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_container(&mut self, container: &mut Container) -> io::Result<usize> {
        read_container(&mut self.inner, container, self.version)
    }

    /// Returns a iterator over records starting from the current stream position.
    ///
    /// The stream is expected to be at the start of a container.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    /// use noodles_fasta as fasta;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// let header = reader.read_header()?;
    ///
    /// for result in reader.records(&header) {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records<'r, 'h: 'r>(&'r mut self, header: &'h sam::Header) -> Records<'r, 'h, R> {
        Records::new(self, header)
    }

    /// Calls `f` for each record, avoiding intermediate `RecordBuf` conversion.
    ///
    /// This is useful when writing CRAM records directly to another format (e.g., BAM), as it
    /// passes each decoded `cram::Record` as `&dyn sam::alignment::Record` without first
    /// converting to `RecordBuf`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// let header = reader.read_header()?;
    ///
    /// reader.for_each_record(&header, |record| {
    ///     // record is &dyn sam::alignment::Record
    ///     Ok(())
    /// })?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn for_each_record<F>(&mut self, header: &sam::Header, mut f: F) -> io::Result<()>
    where
        F: FnMut(&dyn sam::alignment::Record) -> io::Result<()>,
    {
        let mut container = Container::default();

        while self.read_container(&mut container)? != 0 {
            let compression_header = container.compression_header()?;

            for result in container.slices() {
                let slice = result?;
                let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

                let records = slice.records(
                    self.reference_sequence_repository.clone(),
                    header,
                    &compression_header,
                    &core_data_src,
                    &external_data_srcs,
                )?;

                for record in &records {
                    f(record)?;
                }
            }
        }

        Ok(())
    }
}

impl<R> Reader<R>
where
    R: Read + Seek,
{
    /// Seeks the underlying reader to the given position.
    ///
    /// Positions typically come from the associated CRAM index file.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io::{self, SeekFrom};
    /// use noodles_cram as cram;
    /// let mut reader = cram::io::Reader::new(io::empty());
    /// reader.seek(SeekFrom::Start(0))?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.inner.seek(pos)
    }

    /// Returns the current position of the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let mut reader = cram::io::Reader::new(io::empty());
    /// let position = reader.position()?;
    /// assert_eq!(position, 0);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn position(&mut self) -> io::Result<u64> {
        self.inner.stream_position()
    }

    /// Returns an iterator over records that intersects the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram::{self as cram, crai};
    /// use noodles_fasta as fasta;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let header = reader.read_header()?;
    /// let index = crai::fs::read("sample.cram.crai")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<'r, 'h: 'r, 'i: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        index: &'i crai::Index,
        region: &Region,
    ) -> io::Result<Query<'r, 'h, 'i, R>> {
        let reference_sequence_id = header
            .reference_sequences()
            .get_index_of(region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "invalid reference sequence name",
                )
            })?;

        Ok(Query::new(
            self,
            header,
            index,
            reference_sequence_id,
            region.interval(),
        ))
    }

    /// Returns an iterator over unmapped records.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram::{self as cram, crai};
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// let header = reader.read_header()?;
    /// let index = crai::fs::read("sample.cram.crai")?;
    ///
    /// let query = reader.query_unmapped(&header, &index)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn query_unmapped<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        index: &crai::Index,
    ) -> io::Result<impl Iterator<Item = io::Result<sam::alignment::RecordBuf>> + use<'r, 'h, R>>
    {
        let offset = index
            .iter()
            .filter(|record| record.reference_sequence_id().is_none())
            .map(|record| record.offset())
            .min();

        if let Some(offset) = offset {
            self.seek(SeekFrom::Start(offset))?;
        } else {
            self.seek(SeekFrom::End(0))?;
        }

        Ok(self.records(header).filter(|result| {
            result
                .as_ref()
                .map(|record| record.flags().is_unmapped())
                .unwrap_or(true)
        }))
    }
}

impl<R> sam::alignment::io::Read<R> for Reader<R>
where
    R: Read,
{
    fn read_alignment_header(&mut self) -> io::Result<sam::Header> {
        self.read_header()
    }

    fn alignment_records<'a>(
        &'a mut self,
        header: &'a sam::Header,
    ) -> Box<dyn Iterator<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'a> {
        Box::new(
            self.records(header).map(|result| {
                result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
            }),
        )
    }
}
