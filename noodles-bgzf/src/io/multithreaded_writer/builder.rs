use std::{
    io::Write,
    num::NonZero,
    sync::{Arc, atomic::AtomicU64},
};

use bytes::BytesMut;

use super::MultithreadedWriter;
use crate::io::writer::CompressionLevel;

/// A multithreaded BGZF writer builder.
pub struct Builder {
    compression_level: CompressionLevel,
    worker_count: NonZero<usize>,
}

impl Builder {
    /// Sets the compression level.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::{self as bgzf, io::writer::CompressionLevel};
    ///
    /// let builder = bgzf::io::multithreaded_writer::Builder::default()
    ///     .set_compression_level(CompressionLevel::default());
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
    /// use std::num::NonZero;
    /// use noodles_bgzf as bgzf;;
    ///
    /// let builder = bgzf::io::multithreaded_writer::Builder::default()
    ///     .set_worker_count(NonZero::<usize>::MIN);
    /// ```
    pub fn set_worker_count(mut self, worker_count: NonZero<usize>) -> Self {
        self.worker_count = worker_count;
        self
    }

    /// Builds a multithreaded BGZF writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    ///
    /// let builder = bgzf::io::multithreaded_writer::Builder::default()
    ///     .build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> MultithreadedWriter<W>
    where
        W: Write + Send + 'static,
    {
        use super::{State, spawn_deflaters, spawn_writer};

        let worker_count = self.worker_count.get();

        // Position tracking for indexing support
        let position = Arc::new(AtomicU64::new(0));
        let blocks_written = Arc::new(AtomicU64::new(0));

        // Block info channel (unbounded - writer shouldn't block)
        let (block_info_tx, block_info_rx) = crossbeam_channel::unbounded();

        let (write_tx, write_rx) = crossbeam_channel::bounded(worker_count);
        let (deflate_tx, deflate_rx) = crossbeam_channel::bounded(worker_count);

        let writer_handle = spawn_writer(
            writer,
            write_rx,
            Arc::clone(&position),
            Arc::clone(&blocks_written),
            block_info_tx,
        );
        let deflater_handles =
            spawn_deflaters(self.compression_level, self.worker_count, deflate_rx);

        MultithreadedWriter {
            state: State::Running {
                writer_handle,
                deflater_handles,
                write_tx,
                deflate_tx,
                block_info_rx,
            },
            buf: BytesMut::new(),
            current_block_number: 0,
            position,
            blocks_written,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            compression_level: CompressionLevel::default(),
            worker_count: NonZero::<usize>::MIN,
        }
    }
}
