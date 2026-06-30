use std::vec;

use noodles_sam as sam;
use tokio::io::{self, AsyncRead};

use super::Reader;
use crate::io::reader::Container;

/// An async buffered CRAM reader.
pub struct BufReader<R> {
    inner: Reader<R>,
    container: Container,
    records: vec::IntoIter<io::Result<sam::alignment::RecordBuf>>,
}

impl<R> BufReader<R> {
    /// Returns a reference to the underlying reader.
    pub fn get_ref(&self) -> &Reader<R> {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    pub fn get_mut(&mut self) -> &mut Reader<R> {
        &mut self.inner
    }

    /// Returns the underlying reader.
    pub fn into_inner(self) -> Reader<R> {
        self.inner
    }
}

impl<R> BufReader<R>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async buffered CRAM reader.
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
    pub async fn read_record_buf(
        &mut self,
        header: &sam::Header,
        record: &mut sam::alignment::RecordBuf,
    ) -> io::Result<usize> {
        loop {
            if let Some(result) = self.records.next() {
                *record = result?;
                return Ok(1);
            } else if self.fill_buf(header).await? == 0 {
                return Ok(0);
            }
        }
    }

    async fn fill_buf(&mut self, header: &sam::Header) -> io::Result<usize> {
        let container_size = self.inner.read_container(&mut self.container).await?;

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
