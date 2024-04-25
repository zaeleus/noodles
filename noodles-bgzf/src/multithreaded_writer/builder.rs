use std::{io::Write, num::NonZeroUsize};

use super::MultithreadedWriter;
use crate::writer::CompressionLevel;

/// A multithreaded BGZF writer builder.
pub struct Builder {
    compression_level: CompressionLevel,
    worker_count: NonZeroUsize,
}

impl Builder {
    /// Sets the compression level.
    pub fn set_compression_level(mut self, compression_level: CompressionLevel) -> Self {
        self.compression_level = compression_level;
        self
    }

    /// Sets the worker count.
    pub fn set_worker_count(mut self, worker_count: NonZeroUsize) -> Self {
        self.worker_count = worker_count;
        self
    }

    /// Builds a multithreaded BGZF writer from a writer.
    pub fn build_from_writer<W>(self, writer: W) -> MultithreadedWriter<W>
    where
        W: Write + Send + 'static,
    {
        MultithreadedWriter::new(self.compression_level, self.worker_count, writer)
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
