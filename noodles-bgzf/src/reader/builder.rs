use std::io::Read;

use super::{block, Reader};
use crate::Block;

/// A BGZF reader builder.
#[derive(Debug, Default)]
pub struct Builder {
    worker_count: Option<usize>,
}

impl Builder {
    /// Sets the worker count.
    ///
    /// By default, the worker count is set to the number of available logical CPUs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let builder = bgzf::reader::Builder::default().set_worker_count(1);
    /// ```
    pub fn set_worker_count(mut self, worker_count: usize) -> Self {
        self.worker_count = Some(worker_count);
        self
    }

    /// Builds a BGZF reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let reader = bgzf::reader::Builder::default().build_from_reader(io::empty());
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> Reader<R>
    where
        R: Read,
    {
        let worker_count = self.worker_count.unwrap_or_else(num_cpus::get);
        let block_reader = block::Reader::with_worker_count(worker_count, reader);

        Reader {
            inner: block_reader,
            position: 0,
            block: Block::default(),
        }
    }
}
