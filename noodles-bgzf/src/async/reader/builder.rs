use std::num::NonZeroUsize;

use futures::TryStreamExt;
use tokio::io::AsyncRead;

use super::{Inflater, Reader};
use crate::Block;

/// An async BGZF reader builder.
pub struct Builder<R> {
    inner: R,
    worker_count: Option<NonZeroUsize>,
}

impl<R> Builder<R>
where
    R: AsyncRead,
{
    pub(crate) fn new(inner: R) -> Self {
        Self {
            inner,
            worker_count: None,
        }
    }

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
    /// let data = [];
    /// let worker_count = NonZeroUsize::try_from(1)?;
    /// let builder = bgzf::AsyncReader::builder(&data[..]).set_worker_count(worker_count);
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn set_worker_count(mut self, worker_count: NonZeroUsize) -> Self {
        self.worker_count = Some(worker_count);
        self
    }

    /// Builds an async BGZF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let data = [];
    /// let reader = bgzf::AsyncReader::builder(&data[..]).build();
    /// ```
    pub fn build(self) -> Reader<R> {
        let worker_count = self.worker_count.unwrap_or_else(|| {
            // SAFETY: `num_cpus::get` is guaranteed to be non-zero.
            NonZeroUsize::new(num_cpus::get()).unwrap()
        });

        Reader {
            stream: Some(Inflater::new(self.inner).try_buffered(worker_count.get())),
            block: Block::default(),
            position: 0,
            worker_count,
        }
    }
}
