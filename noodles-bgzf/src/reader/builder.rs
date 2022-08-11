use std::{io::Read, num::NonZeroUsize};

use super::{block, Reader};
use crate::Block;

/// A BGZF reader builder.
#[derive(Debug, Default)]
pub struct Builder {
    worker_count: Option<NonZeroUsize>,
}

impl Builder {
    /// Sets the worker count.
    ///
    /// By default, the worker count is set to the number of available logical CPUs.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_bgzf as bgzf;
    ///
    /// let worker_count = NonZeroUsize::try_from(1)?;
    /// let builder = bgzf::reader::Builder::default().set_worker_count(worker_count);
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn set_worker_count(mut self, worker_count: NonZeroUsize) -> Self {
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
        let worker_count = self.worker_count.unwrap_or_else(|| {
            // SAFETY: `num_cpus::get` is guaranteed to be non-zero.
            NonZeroUsize::new(num_cpus::get()).unwrap()
        });

        let block_reader = block::Reader::with_worker_count(worker_count, reader);

        Reader {
            inner: block_reader,
            position: 0,
            block: Block::default(),
        }
    }
}
