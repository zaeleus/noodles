use std::{io::Read, num::NonZeroUsize};

use super::{block, Reader};
use crate::Block;

const DEFAULT_WORKER_COUNT: NonZeroUsize = match NonZeroUsize::new(1) {
    Some(worker_count) => worker_count,
    None => unreachable!(),
};

/// A BGZF reader builder.
#[derive(Debug)]
pub struct Builder {
    worker_count: NonZeroUsize,
}

impl Builder {
    /// Sets the worker count.
    ///
    /// By default, the worker count is set to 1.
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
        self.worker_count = worker_count;
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
        let block_reader = if self.worker_count.get() == 1 {
            block::Inner::Single(block::single::Reader::new(reader))
        } else {
            block::Inner::Multi(block::multi::Reader::with_worker_count(
                self.worker_count,
                reader,
            ))
        };

        Reader {
            inner: block_reader,
            position: 0,
            block: Block::default(),
            gzi: None,
            uncompressed_position: None,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            worker_count: DEFAULT_WORKER_COUNT,
        }
    }
}
