use std::{io::Write, num::NonZeroUsize};

use bytes::BytesMut;

use super::MultithreadedWriter;
use crate::writer::CompressionLevel;

/// A multithreaded BGZF writer builder.
pub struct Builder {
    compression_level: CompressionLevel,
    worker_count: NonZeroUsize,
}

impl Builder {
    /// Sets the compression level.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::{multithreaded_writer::Builder, writer::CompressionLevel};
    /// let builder = Builder::default().set_compression_level(CompressionLevel::default());
    /// ```
    pub fn set_compression_level(mut self, compression_level: CompressionLevel) -> Self {
        self.compression_level = compression_level;
        self
    }

    /// Sets the worker count.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    /// use noodles_bgzf::multithreaded_writer::Builder;
    /// let builder = Builder::default().set_worker_count(NonZeroUsize::MIN);
    /// ```
    pub fn set_worker_count(mut self, worker_count: NonZeroUsize) -> Self {
        self.worker_count = worker_count;
        self
    }

    /// Builds a multithreaded BGZF writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf::multithreaded_writer::Builder;
    /// let builder = Builder::default().build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> MultithreadedWriter<W>
    where
        W: Write + Send + 'static,
    {
        use super::{spawn_deflaters, spawn_writer, State};

        let worker_count = self.worker_count.get();

        let (write_tx, write_rx) = crossbeam_channel::bounded(worker_count);
        let (deflate_tx, deflate_rx) = crossbeam_channel::bounded(worker_count);

        let writer_handle = spawn_writer(writer, write_rx);
        let deflater_handles =
            spawn_deflaters(self.compression_level, self.worker_count, deflate_rx);

        MultithreadedWriter {
            state: State::Running {
                writer_handle,
                deflater_handles,
                write_tx,
                deflate_tx,
            },
            buf: BytesMut::new(),
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            compression_level: CompressionLevel::default(),
            worker_count: NonZeroUsize::MIN,
        }
    }
}
