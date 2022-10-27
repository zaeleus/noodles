use std::num::NonZeroUsize;

use futures::TryStreamExt;
use tokio::io::AsyncRead;

use super::{Inflater, Reader};
use crate::Block;

/// An async BGZF reader builder.
#[derive(Default)]
pub struct Builder {
    worker_count: Option<NonZeroUsize>,
}

impl Builder {
    /// Sets a worker count.
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
    /// let builder = bgzf::r#async::reader::Builder::default()
    ///     .set_worker_count(worker_count);
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn set_worker_count(mut self, worker_count: NonZeroUsize) -> Self {
        self.worker_count = Some(worker_count);
        self
    }

    /// Builds an async BGZF reader with an async reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use tokio::io;
    /// use noodles_bgzf as bgzf;
    /// let reader = bgzf::r#async::reader::Builder::default()
    ///     .build_with_reader(io::empty());
    /// ```
    pub fn build_with_reader<R>(self, reader: R) -> Reader<R>
    where
        R: AsyncRead,
    {
        let worker_count = self.worker_count.unwrap_or_else(|| {
            // SAFETY: `num_cpus::get` is guaranteed to be non-zero.
            NonZeroUsize::new(num_cpus::get()).unwrap()
        });

        Reader {
            stream: Some(Inflater::new(reader).try_buffered(worker_count.get())),
            block: Block::default(),
            position: 0,
            worker_count,
        }
    }
}
