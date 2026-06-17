use std::{
    io::{self, Read},
    vec,
};

use noodles_sam as sam;

use super::{Reader, reader::Container};

/// A buffered CRAM reader.
pub struct BufReader<R> {
    inner: Reader<R>,
    container: Container,
    records: vec::IntoIter<io::Result<sam::alignment::RecordBuf>>,
}

impl<R> BufReader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let reader = cram::io::BufReader::new(cram::io::Reader::new(io::empty()));
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &Reader<R> {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let mut reader = cram::io::BufReader::new(cram::io::Reader::new(io::empty()));
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut Reader<R> {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let reader = cram::io::BufReader::new(cram::io::Reader::new(io::empty()));
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> Reader<R> {
        self.inner
    }
}

impl<R> BufReader<R>
where
    R: Read,
{
    /// Creates a buffered CRAM reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    ///
    /// let mut reader = File::open("sample.cram")
    ///     .map(cram::io::Reader::new)
    ///     .map(cram::io::BufReader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(inner: Reader<R>) -> Self {
        Self {
            inner,
            container: Container::default(),
            records: Vec::new().into_iter(),
        }
    }

    /// Reads a record into an alignment record buffer.
    ///
    /// If successful, 1 is returned. If 0 is returned, the stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram as cram;
    /// use noodles_sam as sam;
    ///
    /// let mut reader = File::open("sample.cram")
    ///     .map(cram::io::Reader::new)
    ///     .map(cram::io::BufReader::new)?;
    ///
    /// let header = reader.get_mut().read_header()?;
    ///
    /// let mut record = sam::alignment::RecordBuf::default();
    ///
    /// while reader.read_record_buf(&header, &mut record)? != 0 {
    ///     // ...
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_record_buf(
        &mut self,
        header: &sam::Header,
        record: &mut sam::alignment::RecordBuf,
    ) -> io::Result<usize> {
        loop {
            if let Some(result) = self.records.next() {
                *record = result?;
                return Ok(1);
            } else if self.fill_buf(header)? == 0 {
                return Ok(0);
            }
        }
    }

    fn fill_buf(&mut self, header: &sam::Header) -> io::Result<usize> {
        let container_size = self.inner.read_container(&mut self.container)?;

        if container_size == 0 {
            return Ok(0);
        }

        let reference_sequence_repository = self.get_ref().reference_sequence_repository().clone();
        let compression_header = self.container.compression_header()?;

        let mut records = Vec::with_capacity(self.container.header().record_count());

        for result in self.container.slices() {
            let slice = result?;

            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

            let record_bufs = slice
                .records(
                    reference_sequence_repository.clone(),
                    header,
                    &compression_header,
                    &core_data_src,
                    &external_data_srcs,
                )?
                .into_iter()
                .map(|record| {
                    sam::alignment::RecordBuf::try_from_alignment_record(header, &record)
                });

            records.extend(record_bufs);
        }

        self.records = records.into_iter();

        Ok(container_size)
    }
}
